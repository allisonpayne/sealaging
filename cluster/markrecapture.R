library(posterior)
library(tidyverse)
library(rstan)

args <- commandArgs(trailingOnly = TRUE)
if (args[1] == "test") {
  test <- TRUE
} else if (args[1] == "full") {
  test <- FALSE
} else {
  stop(sprintf("command line argument must be 'test' or 'full', not %s", args[1]))
}

sealdat <- read_csv("data/raw/128L pull 2023_12_05.csv") %>% 
  filter(age > 3,
         year < 2023)

# int<lower=0> nind;
# int<lower=0> n_occasions;
# array[nind, n_occasions] int<lower=1, upper=3> y;                   

animalIDs <- sort(unique(sealdat$animalID))
years <- sort(unique(sealdat$year))
n_indiv <- length(animalIDs)
n_occ <- length(years)
birth_years <- sealdat %>% 
  arrange(animalID) %>%
  distinct(animalID, yearborn) %>%
  pull(yearborn) - min(years)

sealobs <- matrix(0, nrow = n_indiv, ncol = n_occ)

for (i in seq(n_indiv)) {
  for (j in seq(n_occ)) {
    animalID <- animalIDs[i]
    year <- years[j]
    obs <- sealdat$observed[sealdat$animalID == animalID & sealdat$year == year]
    #workaround for some seals that have multiple observations in a year
    if (length(obs) > 1)
      obs <- "B"
    if (length(obs) == 0)
      sealobs[i, j] <- 3
    else if (obs == "B")
      sealobs[i, j] <- 1 
    else 
      sealobs[i, j] <- 2
  }
}

rownames(sealobs) <- animalIDs
colnames(sealobs) <- years

# Seal resight map
image(years, seq_along(animalIDs), t(sealobs))

# Test parameters
n_indiv_test <- 100
animalIDs_test <- sample(animalIDs, n_indiv_test)

# Setting up directory to save plots
output_path <- args[2]
dir.create(output_path)

# Basic model -------------------------------------------------------------

# Fit the dang model
if (test) {
  markrecapdat <- list(
    nind = n_indiv_test, 
    n_occasions = n_occ, 
    y = sealobs[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod <- stan("cluster/multistate.stan", 
                       data = markrecapdat, 
                       chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat <- list(
    nind = n_indiv, 
    n_occasions = n_occ, 
    y = sealobs
  )
  markrecapmod <- stan("cluster/multistate.stan", 
                       data = markrecapdat, 
                       chains = 4, iter = 2000, cores = 4)
}

# Look at posterior distribution
markrecapdraws <- as_draws_df(markrecapmod) %>% 
  subset_draws(variable = names(markrecapmod)[1:10])
summarize_draws(markrecapdraws)

# Time varying detection  -------------------------------------------------

# Fit the dang model
if (test) {
  markrecapdat2 <- list(
    nind = n_indiv_test, 
    n_occasions = n_occ, 
    y = sealobs[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod2 <- stan("cluster/multistate2.stan", 
                        data = markrecapdat2, 
                        chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat2 <- list(
    nind = n_indiv, 
    n_occasions = n_occ, 
    y = sealobs
  )
  markrecapmod2 <- stan("cluster/multistate2.stan", 
                        data = markrecapdat2, 
                        chains = 4, iter = 2000, cores = 4)
}

# Look at posterior distribution
markrecapdraws2 <- as_draws_df(markrecapmod2)
# Detection probability for breeders and non-breeders
p1 <- markrecapdraws2 %>% 
  select(starts_with(c("pA[", "pB["))) %>% 
  pivot_longer(everything(), 
               names_to = "year", 
               values_to = "detect_prob") %>% 
  mutate(state = ifelse(substr(year, 2, 2) == "A", "Breeder", "Non-breeder"),
         year = years[as.numeric(str_extract(year, "[0-9]+"))]) %>% 
  group_by(year, state) %>% 
  summarize(across(detect_prob, list(mean = mean, 
                                     q025 = \(x) quantile(x, 0.025),
                                     q975 = \(x) quantile(x, 0.975))),
            .groups = "drop") %>% 
  ggplot(aes(year)) +
  geom_ribbon(aes(ymin = detect_prob_q025, ymax = detect_prob_q975),
              fill = "grey80") +
  geom_line(aes(y = detect_prob_mean, color = state), size = 1.5) +
  facet_grid(cols = vars(state)) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c(Breeder = "firebrick", 
                                `Non-breeder` = "cornflowerblue")) +
  labs(x = "Year", y = "Detection probability") +
  theme_classic()
ggsave(file.path(output_path, "DetectProbTimeVary.png"), p1)


# Survival for breeders and non-breeders
p2 <- markrecapdraws2 %>% 
  select(Breeder = phiA, `Non-breeder` = phiB) %>% 
  pivot_longer(everything(),
               names_to = "State",
               values_to = "Survival") %>% 
  ggplot(aes(Survival, color = State)) +
  geom_density(size = 2) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = c(Breeder = "firebrick", 
                                `Non-breeder` = "cornflowerblue")) +
  theme_classic()
ggsave(file.path(output_path, "SurvivalTimeVary.png"), p2)

# Add fixed effect of year on survival and repro --------------------------

# Fit the dang model
if (test) {
  markrecapdat3 <- list(
    nind = n_indiv_test, 
    n_occasions = n_occ, 
    y = sealobs[animalIDs %in% animalIDs_test, ],
    birth_year = birth_years[animalIDs %in% animalIDs_test]
  )
  markrecapmod3 <- stan("cluster/multistate3.stan", 
                        data = markrecapdat3, 
                        chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat3 <- list(
    nind = n_indiv, 
    n_occasions = n_occ, 
    y = sealobs,
    birth_year = birth_years
  )
  markrecapmod3 <- stan("cluster/multistate3.stan", 
                        data = markrecapdat, 
                        chains = 4, iter = 2000, cores = 4)
}

# Look at posterior distribution
markrecapdraws3 <- as_draws_df(markrecapmod3)

# Detection probability for breeders and non-breeders
p3 <- markrecapdraws3 %>% 
  select(starts_with(c("pA[", "pB["))) %>% 
  pivot_longer(everything(), 
               names_to = "year", 
               values_to = "detect_prob") %>% 
  mutate(state = ifelse(substr(year, 2, 2) == "A", "Breeder", "Non-breeder"),
         year = years[as.numeric(str_extract(year, "[0-9]+"))]) %>%
  group_by(year, state) %>% 
  summarize(across(detect_prob, list(mean = mean, 
                                     q025 = \(x) quantile(x, 0.025),
                                     q975 = \(x) quantile(x, 0.975))),
            .groups = "drop") %>% 
  ggplot(aes(year)) +
  geom_ribbon(aes(ymin = detect_prob_q025, ymax = detect_prob_q975),
              fill = "grey80") +
  geom_line(aes(y = detect_prob_mean, color = state), size = 1.5) +
  facet_grid(cols = vars(state)) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c(Breeder = "firebrick", 
                                `Non-breeder` = "cornflowerblue")) +
  labs(x = "Year", y = "Detection probability") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(output_path, "DetectProbFixedEffect.png"), p3)

# Survival by age
inv_logit <- function(x) exp(x) / (1 + exp(x))
p4 <- markrecapdraws3 %>% 
  select(alpha_phiA, alpha_phiB, beta_phiA, beta_phiB) %>% 
  cross_join(tibble(age = 3:22)) %>% 
  mutate(Breeder = inv_logit(alpha_phiA + beta_phiA * age),
         `Non-breeder` = inv_logit(alpha_phiB + beta_phiB * age)) %>% 
  pivot_longer(c(Breeder, `Non-breeder`), 
               names_to = "state", values_to = "phi") %>% 
  group_by(age, state) %>% 
  drop_na(phi) %>% 
  summarize(across(phi, list(mean = mean, 
                             q025 = \(x) quantile(x, 0.025),
                             q975 = \(x) quantile(x, 0.975))),
            .groups = "drop") %>% 
  ggplot(aes(age)) +
  geom_ribbon(aes(ymin = phi_q025, ymax = phi_q975), fill = "grey80") +
  geom_line(aes(y = phi_mean, color = state), size = 1.5) +
  facet_grid(cols = vars(state)) +
  scale_color_manual(values = c(Breeder = "firebrick", 
                                `Non-breeder` = "cornflowerblue")) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(file.path(output_path, "SurvivalFixedEffect.png"), p4)
