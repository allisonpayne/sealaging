library(posterior)
library(rstan)
library(tidyverse)
set.seed(1234)

bern.stan <-
"
data {
  int<lower=0> N;               // number of trials
  int<lower=0, upper=1> y[N];   // success on trial n
}

parameters {
  real<lower=0, upper=1> theta; // chance of success
}

model {
  theta ~ uniform(0, 1);        // prior
  y ~ bernoulli(theta);         // likelihood
}
"

theta <- 0.30
N <- 20
y <- rbinom(N, 1, theta)
print(y)

fit <- stan(model_code = bern.stan,
          data = list(y = y, N = N),
            iter = 5000)

print(fit, probs = c(0.1, 0.9))

theta_draws <- as_draws_df(fit)
plotpostre <- ggplot(theta_draws, aes(theta)) +
geom_density(color = "cornflowerblue", linewidth = 1.5) +
    theme_classic()
ggsave("test.png", plotpostre)

print("DONE")