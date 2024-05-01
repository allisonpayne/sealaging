library(tidyverse)

sealdat <- read_csv(here::here("data/raw/128L pull 2023_12_05.csv")) %>% 
  filter(age > 3)

# int<lower=0> nind;
# int<lower=0> n_occasions;
# array[nind, n_occasions] int<lower=1, upper=3> y;                   

animalIDs <- sort(unique(sealdat$animalID))
years <- sort(unique(sealdat$year))
n_indiv <- length(animalIDs)
n_occ <- length(years)

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

# Seal resight map
image(years, seq_along(animalIDs), t(sealobs))

# Fit the dang model
markrecapdat <- list(
  nind = n_indiv, 
  n_occasions = n_occ, 
  y = sealobs
)
markrecapmod <- rstan::stan("scratch/multistate.stan", 
                            data = markrecapdat, 
                            chains = 1, iter = 200)

# Look at posterior distribution
library(posterior)
markrecapdraws <- as_draws_df(markrecapmod) %>% 
  subset_draws(variable = names(markrecapmod)[1:10])
summarize_draws(markrecapdraws)
