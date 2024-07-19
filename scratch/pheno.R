library(brms)
library(tidyverse)

age_senesce <- 11

sealdat <- read_csv(here::here("data/raw/128L pull 2023_12_05.csv"), 
                    show_col_types = FALSE) %>% 
  mutate(observed = if_else(observed == "B", "Breeder", "Non-breeder"), 
         observed_int = if_else(observed == "Breeder", 1, 0), 
         pup_survived = pupseeneveragain > 0) %>% 
  filter(age > 3, 
         year > 1987) %>% 
  mutate(animalID = factor(animalID),
         year_fct = factor(year),
         age10 = (age - age_senesce) / 10,
         age_cat = factor(age >= age_senesce, 
                          labels = c("Young", "Old"))) %>% 
  group_by(animalID) %>% 
  mutate(firstrepro = ifelse(max(year) < 2020 && any(observed_int == 1),
                             min(age[observed_int == 1]),
                             NA),
         lastrepro = ifelse(max(year) < 2020 && any(observed_int == 1),
                            max(age[observed_int == 1]),
                            NA),
         lastobs = ifelse(max(year) < 2020, max(age), NA)) %>% 
  ungroup()

# Find seals where we have the annual cycle
annual_cycle <- sealdat %>% 
  filter(year >= 2011,
         moltdur >= 7) %>% 
  drop_na(breeddur, tripdur, moltdur) %>% 
  transmute(animalID, year, observed, 
            age, age10, age_cat,
            breeddur, tripdur, moltdur, 
            trip2dur = 365 - breeddur - tripdur - moltdur) %>% 
  mutate(across(ends_with("dur"), list(frac = \(x) x / 365))) %>% 
  distinct()
# Create a matrix column with the fractions
annual_cycle$y <- as.matrix(select(annual_cycle, ends_with("frac")))

pheno_form <- bf(y ~ age10 : age_cat + (1 | animalID) + (1 | year))
durs_typical <- c(breed = 28, postbreed = 74, molt = 37, postmolt = 226)
log_typical <- log(durs_typical) - log(durs_typical[1])
get_prior(pheno_form, data = annual_cycle, family = dirichlet())
pheno_prior <- c(prior_string(str_glue("normal({log_typical['molt']}, 2)"),
                              class = "Intercept",
                              dpar = "mumoltdurfrac"),
                 prior_string(str_glue("normal({log_typical['postbreed']}, 2)"),
                              class = "Intercept",
                              dpar = "mutrip2durfrac"),
                 prior_string(str_glue("normal({log_typical['postmolt']}, 2)"),
                              class = "Intercept",
                              dpar = "mutripdurfrac"))
pheno_fit <- brm(pheno_form, 
                 data = annual_cycle, 
                 family = dirichlet(), 
                 # prior = pheno_prior,
                 iter = 6000, 
                 cores = 4,
                 control = list(adapt_delta = 0.9))
summary(pheno_fit)

pheno_draws <- as_draws_df(pheno_fit, variable = c("b_*"), regex = TRUE)
pheno_grid <- expand_grid(
  age10 = seq(-0.7, 1.0, by = 0.1),
  age_cat = factor(c("Young", "Old"), levels = c("Young", "Old"))
) %>% 
  mutate(age = age10 * 10 + age_senesce) %>% 
  filter((age_cat == "Young" & age < age_senesce) |
           (age_cat == "Old" & age >= age_senesce))

softmax <- \(xi, x) exp(xi) / rowSums(exp(x))
eta <- function(age_cat, age10, alpha, beta_old, beta_young) {
  alpha + ifelse(age_cat == "Old", beta_old, beta_young) * age10
}
pheno_pred <- pheno_draws %>% 
  cross_join(pheno_grid) %>% 
  mutate(
    eta_breed = 0, # reference level
    eta_trip = eta(age_cat, age10, 
                   b_mutripdurfrac_Intercept,
                   `b_mutripdurfrac_age10:age_catOld`,
                   `b_mutripdurfrac_age10:age_catYoung`),
    eta_molt = eta(age_cat, age10, 
                   b_mumoltdurfrac_Intercept,
                   `b_mumoltdurfrac_age10:age_catOld`,
                   `b_mumoltdurfrac_age10:age_catYoung`),
    eta_trip2 = eta(age_cat, age10, 
                    b_mutrip2durfrac_Intercept,
                    `b_mutrip2durfrac_age10:age_catOld`,
                    `b_mutrip2durfrac_age10:age_catYoung`),
    breed = softmax(eta_breed, cbind(eta_breed, eta_trip, eta_molt, eta_trip2)),
    trip = softmax(eta_trip, cbind(eta_breed, eta_trip, eta_molt, eta_trip2)),
    molt = softmax(eta_molt, cbind(eta_breed, eta_trip, eta_molt, eta_trip2)),
    trip2 = softmax(eta_trip2, cbind(eta_breed, eta_trip, eta_molt, eta_trip2)),
    across(breed:trip2, \(x) x * 365)
) %>% 
  as_tibble()
  
# summarized raw data
pheno_summ <- annual_cycle %>% 
  select(breed = breeddur,
         trip = tripdur,
         molt = moltdur,
         trip2 = trip2dur,
         age, age_cat) %>% 
  pivot_longer(breed:trip2, names_to = "phase", values_to = "dur") %>% 
  group_by(age, age_cat, phase) %>% 
  summarize(across(dur, list(mean = mean, se = \(x) sd(x) / sqrt(length(x)))),
            .groups = "drop")

pheno_pred %>% 
  pivot_longer(breed:trip2, names_to = "phase", values_to = "dur") %>% 
  group_by(age, age_cat, phase) %>% 
  summarize(across(dur, list(mean = mean, 
                             p025 = \(x) quantile(x, 0.025), 
                             p975 = \(x) quantile(x, 0.975))),
            .groups = "drop") %>% 
  ggplot(aes(age)) +
  # model
  geom_ribbon(aes(fill = age_cat,
                  ymin = dur_p025,
                  ymax = dur_p975),
              alpha = 0.3) +
  geom_line(aes(color = age_cat, y = dur_mean)) +
  # raw
  geom_pointrange(aes(y = dur_mean, 
                      ymin = dur_mean - 1.96 * dur_se, 
                      ymax = dur_mean + 1.96 * dur_se),
                  pheno_summ) +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_color_manual(values = c("#7fbc41", "#de77ae")) +
  scale_fill_manual(values = c("#7fbc41", "#de77ae")) +
  labs(x = "Age", y = "Duration (days)") +
  facet_wrap(~phase, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none")

