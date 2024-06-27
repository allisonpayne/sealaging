library(lme4)
library(bayesplot)
library(patchwork)
library(ggeffects)
library(posterior)
library(rstan)
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
  mutate(longevity = ifelse(max(year) < 2020,
                            max(age),
                            NA),
         longevity10 = longevity / 10,
         recruit_age = ifelse(max(year) < 2020 && any(observed_int == 1),
                              min(age[observed_int == 1]),
                              NA)) %>% 
  ungroup()

sealdat2 <- sealdat %>%
  select(animalID, year, age, observed_int) %>%
  mutate(across(c(animalID, year), \(x) as.integer(factor(x))))
fit <- stan("scratch/sealaging.stan", "sealaging",
            data = list(N = nrow(sealdat2),
                        A = n_distinct(sealdat2$animalID),
                        Y = n_distinct(sealdat2$year),
                        b = sealdat2$observed_int,
                        k = sealdat2$age,
                        a = sealdat2$animalID,
                        y = sealdat2$year),
            control = list(adapt_delta = 0.8),
            iter = 10000, chains = 4, cores = 4)

# Diagnostics
print(fit, pars = c("alpha1", "beta_k", "tau1", "tau2", "sigma_u", "sigma_v"))

# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
draws_fit <- as_draws_df(fit)
draws_fit2 <- subset_draws(
  draws_fit, 
  variable = c("alpha1", "beta_k", "tau1", "tau2", "sigma_u", "sigma_v")
)
lp_fit <- log_posterior(fit)
np_fit <- nuts_params(fit)

# 0/40000 divergent transitions (0%) Woo!
np_fit %>% 
  filter(Parameter == "divergent__") %>% 
  summarize(total = sum(Value), proportion = mean(Value))

color_scheme_set("darkgray")
mcmc_trace(draws_fit,
           pars = c("tau1", "tau2", "beta_k[1]", "beta_k[2]", "beta_k[3]"))
mcmc_parcoord(draws_fit, np = np_fit,
              pars = c("tau1", "tau2", "beta_k[1]", "beta_k[2]", "beta_k[3]"))
mcmc_pairs(draws_fit, np = np_fit, 
           pars = c("tau1", "tau2", "beta_k[1]", "beta_k[2]", "beta_k[3]"),
           off_diag_args = list(size = 0.75))

inv_logit <- \(x) exp(x) / (1 + exp(x))

## Where are the taus?
tau_corr <- draws_fit2 %>% 
  pivot_longer(c(tau1, tau2), names_to = "breakpoint", values_to = "age")
ggplot(tau_corr, aes(age)) +
  geom_density(aes(color = breakpoint, fill = breakpoint), alpha = 0.5) +
  # geom_vline(xintercept = mean(foo$tau), color = "firebrick") +
  theme_classic()

breed_summ <- sealdat %>% 
  group_by(age) %>% 
  summarize(p_mean = mean(observed_int),
            p_se = sqrt(p_mean * (1 - p_mean) / n()),
            p_lwr = p_mean - p_se,
            p_upr = p_mean + p_se)
age_eta <- function(age, a, b1, b2, b3, t1, t2) {
  beta <- case_when(
    age < t1 ~ b1,
    age < t2 ~ b2,
    TRUE ~ b3
  )
  a + beta * age
}
cross_join(draws_fit2, tibble(age = 4:22)) %>% 
  mutate(age0 = (age - mean(sealdat$age)) / sd(sealdat$age),
         eta = age_eta(age0, alpha1, 
                       `beta_k[1]`, `beta_k[2]`, `beta_k[3]`, 
                       tau1, tau2),
         p = inv_logit(eta)) %>% 
  group_by(age) %>% 
  summarize(p_mean = mean(p),
            p_025 = quantile(p, 0.025),
            p_975 = quantile(p, 0.975)) %>% 
  ggplot(aes(age)) +
  geom_ribbon(aes(ymin = p_025, ymax = p_975), alpha = 0.25) +
  geom_line(aes(y = p_mean)) +
  geom_pointrange(aes(y = p_mean, ymin = p_lwr, ymax = p_upr),
                  breed_summ) + 
  theme_classic()

## Are we missing animals more or less with age?
missing_seals <- sealdat %>% 
  # Need to see a seal 3+ years for calculation to make sense
  group_by(animalID) %>% 
  mutate(n_years = n_distinct(year)) %>% 
  ungroup() %>% 
  filter(n_years >= 3,
         age <= 20) %>% 
  # Fill in non-resighted years
  group_by(animalID) %>% 
  arrange(year) %>% 
  reframe(tibble(year2 = (min(year) + 1):(max(year) - 1),
                 sighted = as.integer(year2 %in% year),
                 age = year2 - yearborn[1],
                 age10 = (age - age_senesce) / 10)) %>% 
  ungroup() %>% 
  rename(year = year2) %>% 
  mutate(age_cat = factor(age >= age_senesce, 
                          labels = c("Young", "Old")),
         year_fct = factor(year)) %>% 
  left_join(select(sealdat, animalID, year, observed_int),
            by = c("animalID", "year"))

missing_mod <- lme4::glmer(sighted ~ age + (1 | year), 
                           missing_seals, 
                           family = binomial)
summary(missing_mod)
missing_pred <- ggeffects::ggpredict(missing_mod, 
                                     terms = "age [all]") %>% 
  as_tibble() %>% 
  rename(age = x,
         sighted = predicted)

missing_seals %>% 
  group_by(age) %>% 
  summarize(sighted_mean = mean(sighted),
            sighted_se = sqrt(sighted_mean * (1 - sighted_mean) / n()),
            sighted_lwr = sighted_mean - sighted_se,
            sighted_upr = sighted_mean + sighted_se) %>% 
  ggplot(aes(age)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              missing_pred,
              alpha = 0.2) +
  geom_line(aes(y = sighted), 
            missing_pred, 
            color = "cornflowerblue",
            linewidth = 2) +
  geom_pointrange(aes(y = sighted_mean, 
                      ymin = sighted_lwr, ymax = sighted_upr)) +
  scale_x_continuous("Age") +
  scale_y_continuous("Resight probability", 
                     labels = scales::label_percent(accuracy = 1)) +
  theme_classic()

## Are the missing seals affecting our inferences?
# If unobserved seals are less likely to be breeders, how big does the effect
# size have to be to change our results?

h1_mod <- glmer(
  observed_int ~ age10 : age_cat + (1 | animalID) + (1 | year_fct),
  sealdat,
  family = "binomial",
  control = glmerControl(optimizer = "bobyqa")
)

missing_seals2 <- sealdat %>% 
  group_by(animalID) %>% 
  mutate(n_years = n_distinct(year)) %>% 
  ungroup() %>% 
  # Fill in non-resighted years
  group_by(animalID) %>% 
  arrange(year) %>% 
  reframe(tibble(year2 = min(year):max(year),
                 sighted = as.integer(year2 %in% year),
                 age = year2 - yearborn[1],
                 age10 = (age - age_senesce) / 10)) %>% 
  ungroup() %>% 
  rename(year = year2) %>% 
  mutate(age_cat = factor(age >= age_senesce, 
                          labels = c("Young", "Old")),
         year_fct = factor(year)) %>% 
  left_join(select(sealdat, animalID, year, observed_int),
            by = c("animalID", "year")) %>% 
  mutate(p_h1 = predict(h1_mod, newdata = ., type = "response"))

impute_missing <- function(dat, offset) {
  # Subtract eff_sz * se(intercept) from the model intercept 
  h1a_mod <- h1_mod
  h1a_mod@beta[1] <- h1a_mod@beta[1] + offset
  
  # Impute missing values
  impute <- \(obs, p) coalesce(obs, rbinom(length(p), 1, p))
  dat_impute <- dat %>% 
    mutate(p_h1a = predict(h1a_mod, newdata = ., type = "response"),
           # Impute missing values
           obs_h1 = impute(observed_int, p_h1),
           obs_h1a = impute(observed_int, p_h1a))
  
  # Fit senescence models
  fit_mod <- function(y) {
    # Replace observed_int with imputed observation
    glmer(
      update(observed_int ~ age10 : age_cat + (1 | animalID) + (1 | year_fct),
             as.formula(paste(y, "~."))),
      dat_impute,
      family = "binomial",
      control = glmerControl(optimizer = "bobyqa")
    )
  }
  list(
    h1_mod = fit_mod("obs_h1"),
    h1a_mod = fit_mod("obs_h1a")
  )
}

h1_sens <- tibble(
  offset = rep(seq(0, -6, by = -0.5), each = 100),
  iter = seq(length(offset)),
  result = map(offset, \(o) impute_missing(missing_seals2, o))
)

h1_sens_df <- h1_sens %>% 
  rowwise() %>% 
  reframe({
    h1_result <- summary(result$h1_mod)$coefficients %>% 
      as_tibble(rownames = "Param") %>% 
      pivot_longer(-Param, names_to = "metric") %>%
      mutate(model = "h1")
    h1a_result <- summary(result$h1a_mod)$coefficients %>% 
      as_tibble(rownames = "Param") %>% 
      pivot_longer(-Param, names_to = "metric") %>%
      mutate(model = "h1a")
    rbind(h1_result, h1a_result) %>% 
      mutate(offset = offset, iter = iter)
  }) %>% 
  relocate(iter, offset)

h1_sens_df %>% 
  filter(Param %in% c("age10:age_catYoung", "age10:age_catOld"),
         metric %in% c("Estimate", "Pr(>|z|)")) %>% 
  pivot_wider(names_from = "metric", values_from = "value") %>% 
  ggplot(aes(Estimate, color = interaction(model, Param))) +
  geom_boxplot() +
  facet_wrap(~ offset) +
  theme_classic()

impute_df <- h1_sens %>% 
  group_by(offset, iter) %>% 
  reframe({
    h1 <- result[[1]]$h1_mod@frame %>% 
      mutate(model = "h1") %>% 
      rename(imputed_obs = obs_h1)
    h1a <- result[[1]]$h1a_mod@frame %>% 
      mutate(model = "h1a") %>% 
      rename(imputed_obs = obs_h1a)
    rbind(h1, h1a)
  }) %>% 
  mutate(age = age10 * 10 + age_senesce)

impute_mods <- impute_df %>% 
  filter(model == "h1a") %>% 
  group_by(offset) %>% 
  group_map(\(x, y) {
    glmer(
      imputed_obs ~ age10 : age_cat + (1 | animalID) + (1 | year_fct),
      x,
      family = "binomial",
      control = glmerControl(optimizer = "bobyqa")
    )
  })

impute_summ <- impute_df %>% 
  filter(model == "h1a") %>% 
  group_by(offset, age) %>% 
  summarize(n_observed = n(), 
            perc_breed = sum(imputed_obs) / n_observed, 
            se = (perc_breed * (1 - perc_breed) / n_observed)^0.5,
            .groups = "drop")

ggplot(filter(impute_summ,, aes(age, perc_breed, color = factor(offset))) +
  geom_point() +
  geom_smooth(aes(color = factor(offset)), se = FALSE) +
  facet_wrap(~model) +
  theme_classic()

example <- h1_sens %>% 
  filter(offset == -6) %>% 
  pull(result) %>% 
  first()

h1_dat <- example$h1_mod@frame
h1a_dat <- example$h1a_mod@frame

h1_dat_summ <- h1_dat %>% 
  mutate(age = age10 * 10 + age_senesce) %>% 
  group_by(age) %>% 
  summarize(n_observed = n(), 
            perc_breed = sum(obs_h1) / n_observed, 
            se = (perc_breed * (1 - perc_breed) / n_observed)^0.5)

h1_predict <- ggeffects::ggpredict(example$h1_mod, 
                                   terms = c("age10 [all]", "age_cat")) %>% 
  as_tibble() %>% 
  mutate(age = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old"))) %>% 
  filter((age_cat == "Young" & age < age_senesce) |
           (age_cat == "Old" & age >= age_senesce),
         age <= 20)

p_h1 <- ggplot(h1_predict, aes(age, predicted)) +
  geom_ribbon(aes(fill = age_cat, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2) +
  geom_line(aes(color = age_cat), linewidth = 1.2) +
  geom_pointrange(aes(y = perc_breed, 
                      ymin = perc_breed - se, 
                      ymax = perc_breed + se),
                  h1_dat_summ) +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous("Breeding", labels = scales::percent) +
  scale_color_manual(values = c("#7fbc41", "#de77ae")) +
  scale_fill_manual(values = c("#7fbc41", "#de77ae")) +
  labs(x = "Age (years)",
       title = "Offset = 0") +
  theme_classic() +
  theme(legend.position = "none") 

h1a_dat_summ <- h1a_dat %>% 
  mutate(age = age10 * 10 + age_senesce) %>% 
  group_by(age) %>% 
  summarize(n_observed = n(), 
            perc_breed = sum(obs_h1a) / n_observed, 
            se = (perc_breed * (1 - perc_breed) / n_observed)^0.5)

h1a_predict <- ggeffects::ggpredict(example$h1a_mod, 
                                    terms = c("age10 [all]", "age_cat")) %>% 
  as_tibble() %>% 
  mutate(age = x * 10 + age_senesce,
         age_cat = factor(group, levels = c("Young", "Old"))) %>% 
  filter((age_cat == "Young" & age < age_senesce) |
           (age_cat == "Old" & age >= age_senesce),
         age <= 20)

p_h1a <- ggplot(h1a_predict, aes(age, predicted)) +
  geom_ribbon(aes(fill = age_cat, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2) +
  geom_line(aes(color = age_cat), linewidth = 1.2) +
  geom_pointrange(aes(y = perc_breed, 
                      ymin = perc_breed - se, 
                      ymax = perc_breed + se),
                  h1a_dat_summ) +
  geom_vline(xintercept = age_senesce - 0.5, linetype = "dashed") +
  scale_y_continuous("Breeding", labels = scales::percent) +
  scale_color_manual(values = c("#7fbc41", "#de77ae")) +
  scale_fill_manual(values = c("#7fbc41", "#de77ae")) +
  labs(x = "Age (years)",
       title = "Offset = -6") +
  theme_classic() +
  theme(legend.position = "none") 

p_h1 | p_h1a

intercept0 <- summary(h1_mod)$coef["(Intercept)", "Estimate"]
h1_sens_df %>% 
  filter(Param %in% c("age10:age_catYoung", "age10:age_catOld")) %>% 
  mutate(age_cat = ifelse(Param == "age10:age_catYoung", "yng", "old")) %>%
  pivot_wider(names_from = metric, values_from = value) %>% 
  group_by(offset, model, iter) %>% 
  summarize(young_nonneg = Estimate[age_cat == "yng"] - 
              1.96 * `Std. Error`[age_cat == "yng"] >= 0,
            old_neg = Estimate[age_cat == "old"] + 
              1.96 * `Std. Error`[age_cat == "old"] < 0,
            .groups = "drop_last") %>% 
  summarize(h1_evidence = mean(young_nonneg & old_neg),
            young_nonneg = mean(young_nonneg),
            old_neg = mean(old_neg),
            .groups = "drop") %>% 
  mutate(intercept = inv_logit(intercept0 + offset)) %>% 
  pivot_longer(c(young_nonneg, old_neg, h1_evidence),
               names_to = "evidence",
               values_to = "proportion") %>% 
  ggplot(aes(intercept, proportion, color = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~evidence) +
  theme_bw() +
  theme(strip.placement = "outside")
