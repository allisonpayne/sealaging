library(brms)
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
  ungroup() %>% 
  mutate(age2 = (age - 4) / 10) 

bform <- bf(
  # intercept and fixed effect of longevity
  observed_int ~ b0 + bl * longevity +
    # pre-change slope
    b1 * (age2 - omega) * step(omega - age2) +
    # post-change slope
    b2 * (age2 - omega) * step(age2 - omega),
  # intercept varies by animal and year
  b0 ~ 1 + (1 | animalID + year),
  # estimate longevity coef and change point for population
  bl + b1 + b2 + omega ~ 1,
  nl = TRUE
)

bprior <- prior(normal(0, 3), nlpar = "b0") +
  prior(normal(0, 3), nlpar = "b1") +
  prior(normal(0, 3), nlpar = "b2") +
  prior(normal(0, 3), nlpar = "bl") + 
  prior(uniform(0, 1.8), nlpar = "omega")

fit <- brm(bform, data = sealdat, prior = bprior, chains = 4, iter = 10000)
