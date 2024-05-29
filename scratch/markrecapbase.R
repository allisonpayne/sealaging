library(posterior)
library(tidyverse)
library(rstan)

inv_logit <- \(x) 1 / (1 + exp(-x))

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

# Read command line
args <- commandArgs(trailingOnly = TRUE)
if (args[1] == "test") {
  test <- TRUE
} else if (args[1] == "full") {
  test <- FALSE
} else {
  stop(sprintf("First argument must be 'test' or 'full', not %s", args[1]))
}

# Setting up directory to save plots
output_path <- args[2]
dir.create(output_path)

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
  markrecapmod0 <- stan("scratch/sealaging.stan",
                       data = markrecapdat,
                       chains = 4, iter = 8000, cores = 4)
}

markrecapdraws <- as_draws_df(markrecapmod)
summarize_draws(markrecapdraws)
stan_trace(markrecapmod, 
           pars = c("pb_maturity", "lps[1]", "lps[2]", "lpb[1]", "lpb[2]",
                    "lpd_mu[1]", "lpd_mu[2]"))

# Fit time-varying detection model ----------------------------------------
if (test) {
  markrecapdat <- list(
    n_ind = n_ind_test,
    n_occ = n_occ,
    obs = obs[animalIDs %in% animalIDs_test, ],
    age = age[animalIDs %in% animalIDs_test, ]
  )
  markrecapmod <- stan("scratch/sealaging3.stan",
                       data = markrecapdat,
                       chains = 1, iter = 200, cores = 1)
} else {
  markrecapdat <- list(
    n_ind = n_ind,
    n_occ = n_occ,
    obs = obs,
    age = age
  )
  markrecapmod <- stan("scratch/sealaging3.stan",
                       data = markrecapdat,
                       chains = 4, iter = 8000, warmup = 2000, cores = 4)
}

markrecapdraws <- as_draws_df(markrecapmod)
summarize_draws(markrecapdraws)
stan_trace(markrecapmod, 
           pars = c("logit_ps[1]", "logit_ps[2]", 
                    "logit_pb[1]", "logit_pb[2]", 
                    "logit_pd_mu[1]", "logit_pd_mu[2]"))
# Detection probability
pd_tbl <- markrecapdraws %>% 
  select(matches("^logit_pd\\[.*\\]")) %>% 
  pivot_longer(everything(), names_to = "param", values_to = "logit_pd") %>% 
  mutate(state = factor(as.numeric(str_extract(param, "pd\\[([0-9]+),[0-9]+\\]", 1)),
                        labels = c("Breeding", "Non-breeding")),
         t = as.numeric(str_extract(param, "pd\\[[0-9]+,([0-9]+)\\]", 1)),
         year = years[t],
         pd = inv_logit(logit_pd))
pd_tbl %>% 
  group_by(state, year) %>% 
  summarize(across(pd, list(mean = mean, 
                            q025 = \(x) quantile(x, 0.025), 
                            q975 = \(x) quantile(x, 0.975))),
            .groups = "drop") %>% 
  ggplot(aes(year)) +
  geom_ribbon(aes(ymin = pd_q025, ymax = pd_q975), fill = "grey80") +
  geom_line(aes(y = pd_mean, color = state), linewidth = 1.5) +
  scale_color_manual(values = c("cornflowerblue", "firebrick")) +
  scale_y_continuous("p(Detection)", labels = scales::percent) +
  facet_grid(rows = vars(state)) +
  theme_classic() +
  theme(legend.position = c(0.1, 0.6),
        legend.justification = c(0, 0.5),
        legend.title = element_blank()) 
# Survival
markrecapdraws %>% 
  select(matches("logit_ps")) %>% 
  pivot_longer(everything(), names_to = "param", values_to = "logit_ps") %>% 
  mutate(state = factor(as.numeric(substr(param, 10, 10)),
                        labels = c("Breeding", "Non-breeding")),
         ps = inv_logit(logit_ps)) %>% 
  ggplot(aes(ps)) +
  geom_density(aes(color = state), linewidth = 1.5) +
  scale_color_manual(values = c("cornflowerblue", "firebrick")) +
  scale_x_continuous("p(Survival)", labels = scales::percent) +
  theme_classic() +
  theme(legend.position = c(0.1, 0.9),
        legend.justification = c(0, 1),
        legend.title = element_blank())
