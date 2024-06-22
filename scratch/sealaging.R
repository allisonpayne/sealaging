library(posterior)
library(rstan)
library(tidyverse)

sealdat <- read_csv(here::here("data/raw/128L pull 2023_12_05.csv"), 
                    show_col_types = FALSE) %>% 
  mutate(observed = if_else(observed == "B", "Breeder", "Non-breeder"), 
         observed_int = if_else(observed == "Breeder", 1, 0), 
         pup_survived = pupseeneveragain > 0) %>% 
  filter(age > 3, 
         year > 1987) %>% 
  mutate(animalID = factor(animalID),
         year_fct = factor(year),
         age10 = age / 10) %>% 
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
  select(animalID, year, age, observed_int, longevity) %>% 
  drop_na(longevity) %>% 
  mutate(across(c(age, longevity), \(x) (x - mean(x)) / sd(x)),
         across(c(animalID, year), \(x) as.integer(factor(x))))
fit <- stan("scratch/sealaging.stan", "sealaging",
            data = list(N = nrow(sealdat2),
                        A = n_distinct(sealdat2$animalID),
                        Y = n_distinct(sealdat2$year),
                        b = sealdat2$observed_int,
                        k = sealdat2$age,
                        a = sealdat2$animalID,
                        y = sealdat2$year,
                        l = sealdat2$longevity),
            iter = 8000, chains = 4, cores = 4)

fit_draws <- subset_draws(
  as_draws_df(fit), 
  variable = c("alpha", "beta_k", "beta_l", "tau", "sigma_u", "sigma_v")
)

inv_logit <- \(x) exp(x) / (1 + exp(x))

## Where is tau?
foo <- fit_draws %>% 
  mutate(tau = tau * sd(sealdat$age) + mean(sealdat$age))
ggplot(foo, aes(tau)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = mean(foo$tau), color = "firebrick") +
  theme_classic()

breed_summ <- sealdat %>% 
  group_by(age) %>% 
  summarize(p_mean = mean(observed_int),
            p_se = sqrt(p_mean * (1 - p_mean) / n()),
            p_lwr = p_mean - p_se,
            p_upr = p_mean + p_se)
cross_join(fit_draws, tibble(age = 4:22)) %>% 
  mutate(age0 = (age - mean(sealdat$age)) / sd(sealdat$age),
         eta = alpha + ifelse(age0 < tau, `beta_k[1]`, `beta_k[2]`) * age0,
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

sealdat %>% 
  group_by(age) %>% 
  summarize(breeding = mean(observed_int)) %>% 
  ggplot(aes(age, breeding)) +
  geom_line()
