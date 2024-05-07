library(posterior)
library(tidyverse)
library(rstan)

sealdat <- read_csv("data/raw/128L pull 2023_12_05.csv") %>% 
  filter(age > 3,
         year < 2023,
         yearborn >= min(year) - 2)            

animalIDs <- sort(unique(sealdat$animalID))
years <- sort(unique(sealdat$year))
n_ind <- length(animalIDs)
n_occ <- length(years)

obs <- matrix(0, nrow = n_ind, ncol = n_occ)
age <- matrix(0, nrow = n_ind, ncol = n_occ)

for (i in seq(n_ind)) {
  animalID <- animalIDs[i]
  yearborn <- sealdat$yearborn[sealdat$animalID == animalID][1]
  for (j in seq(n_occ)) {
    year <- years[j]
    obsij <- sealdat$observed[sealdat$animalID == animalID & sealdat$year == year]
    #workaround for some seals that have multiple observations in a year
    if (length(obsij) > 1)
      obsij <- "B"
    if (length(obsij) == 0)
      obs[i, j] <- 3
    else if (obsij == "B")
      obs[i, j] <- 1 
    else 
      obs[i, j] <- 2
    age[i, j] <- year - yearborn
  }
}

rownames(obs) <- rownames(age) <- animalIDs
colnames(obs) <- colnames(age) <- years

# Seal resight map
image(years, seq_along(animalIDs), t(obs))

# Test parameters
n_ind_test <- 100
animalIDs_test <- sample(animalIDs, n_ind_test)

# Setting up directory to save plots
# output_path <- args[2]
# dir.create(output_path)

# Fit model ---------------------------------------------------------------
if (test) {
  markrecapdat <- list(
    n_ind = n_ind_test,
    n_occ = n_occ,
    obs = obs[animalIDs %in% animalIDs_test, ],
    age = age[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod <- stan("scratch/sealaging.stan",
                       data = markrecapdat,
                       chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat <- list(
    n_ind = n_ind,
    n_occ = n_occ,
    obs = obs,
    age = age
  )
  markrecapmod <- stan("scratch/sealaging.stan",
                       data = markrecapdat,
                       chains = 4, iter = 8000, cores = 4)
}

markrecapdraws <- as_draws_df(markrecapmod)
summarize_draws(markrecapdraws) %>% view()
stan_trace(markrecapmod, 
           pars = c("pb_maturity", "ps[1]", "ps[2]", "pd[1]", "pd[2]", 
                    "theta12[1,1]", "theta12[2,1]"))

# Basic model -------------------------------------------------------------

print("BEGINNING BASIC MODEL")

# Fit the dang model
if (test) {
  markrecapdat <- list(
    nind = n_indiv_test, 
    n_occasions = n_occ, 
    y = obs[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod <- stan("cluster/multistate.stan", 
                       data = markrecapdat, 
                       chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat <- list(
    nind = n_indiv, 
    n_occasions = n_occ, 
    y = obs
  )
  markrecapmod <- stan("cluster/multistate.stan", 
                       data = markrecapdat, 
                       chains = 4, iter = 2000, cores = 4)
}

# Look at posterior distribution
markrecapdraws <- as_draws_df(markrecapmod) %>% 
  subset_draws(variable = names(markrecapmod)[1:10])
summarize_draws(markrecapdraws)

print("BASIC MODEL COMPLETE")

# Time varying detection  -------------------------------------------------

print("BEGINNING TIME-VARYING DETECTION MODEL")

# Fit the dang model
if (test) {
  markrecapdat2 <- list(
    nind = n_indiv_test, 
    n_occasions = n_occ, 
    y = obs[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod2 <- stan("cluster/multistate2.stan", 
                        data = markrecapdat2, 
                        chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat2 <- list(
    nind = n_indiv, 
    n_occasions = n_occ, 
    y = obs
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

print("TIME-VARYING DETECTION MODEL COMPLETE")

# Add fixed effect of year on survival and repro --------------------------

print("BEGINNING AGE FIXED EFFECT MODEL")

# Fit the dang model
if (test) {
  markrecapdat3 <- list(
    nind = n_indiv_test, 
    n_occasions = n_occ, 
    y = obs[animalIDs %in% animalIDs_test, ],
    birth_year = birth_years[animalIDs %in% animalIDs_test]
  )
  markrecapmod3 <- stan("cluster/multistate3.stan", 
                        data = markrecapdat3, 
                        chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat3 <- list(
    nind = n_indiv, 
    n_occasions = n_occ, 
    y = obs,
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

print("AGE FIXED EFFECT MODEL COMPLETE")