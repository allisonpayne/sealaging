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

# Fit basic model ---------------------------------------------------------
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
                    "pb[1]", "pb[2]"))

# Fit time-varying detection model ----------------------------------------
if (test) {
  markrecapdat <- list(
    n_ind = n_ind_test,
    n_occ = n_occ,
    obs = obs[animalIDs %in% animalIDs_test, ],
    age = age[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod <- stan("scratch/sealaging2.stan",
                       data = markrecapdat,
                       chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat <- list(
    n_ind = n_ind,
    n_occ = n_occ,
    obs = obs,
    age = age
  )
  markrecapmod <- stan("scratch/sealaging2.stan",
                       data = markrecapdat,
                       chains = 4, iter = 8000, cores = 4)
}

markrecapdraws <- as_draws_df(markrecapmod)
summarize_draws(markrecapdraws) %>% view()
stan_trace(markrecapmod, 
           pars = c("pb_maturity", "ps[1]", "ps[2]", "pd[1,1]", "pd[2,1]", 
                    "pb[1]", "pb[2]"))
pd_tbl <- markrecapdraws %>% 
  select(matches("^pd\\[.*\\]")) %>% 
  pivot_longer(everything(), names_to = "param", values_to = "pd") %>% 
  mutate(state = as.numeric(str_extract(param, "pd\\[([0-9]+),[0-9]+\\]", 1)),
         t = as.numeric(str_extract(param, "pd\\[[0-9]+,([0-9]+)\\]", 1)),
         year = years[t])
pd_tbl %>% 
  group_by(state, year) %>% 
  summarize(across(pd, list(mean = mean, 
                            q025 = \(x) quantile(x, 0.025), 
                            q975 = \(x) quantile(x, 0.975))),
            .groups = "drop") %>% 
  ggplot(aes(year)) +
  geom_ribbon(aes(ymin = pd_q025, ymax = pd_q975), fill = "grey80") +
  geom_line(aes(y = pd_mean), color = "cornflowerblue") +
  facet_grid(rows = vars(state)) +
  theme_classic()
