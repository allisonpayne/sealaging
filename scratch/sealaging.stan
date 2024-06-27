data {
  int N;          // number of obs (seal x year)
  int A;          // number of distinct animals
  int Y;          // number of distinct years

  array[N] int b; // breeding status (0 = non-breeding, 1 = breeding)
  vector[N] k;    // ages
  array[N] int a; // animal IDs
  array[N] int y; // year IDs
  // vector[N] l;    // longevities
}
transformed data {
  real mean_k = mean(k);
  real sd_k = sd(k);
  vector[N] k_std = (k - mean_k) / sd_k;
}
parameters {
  vector[A] u;                 // Ranef of animal
  vector[Y] v;                 // Ranef of year
  real<lower=0> sigma_u;       // SD of animal ranef
  real<lower=0> sigma_v;       // SD of year ranef
  real alpha1;                 // Intercept 1 (i.e. p(breeding) at tau[1])
  vector[3] beta_k;            // Age coef (dev, prime, senesce)
  // real beta_l;           // Coefficient for longevity
  real<lower=5,upper=9> tau1;   // Breakpoint dev -> prime
  real<lower=10,upper=17> tau2; // Breakpoint prime -> senescent
}
model {
  // Standardize tau
  real tau1_std = (tau1 - mean_k) / sd_k;
  real tau2_std = (tau2 - mean_k) / sd_k;
  // Prior on alpha1 is the average breeding probability
  real alpha1_mu = logit(mean(b));
  // Calculate alpha2
  real alpha2 = beta_k[2] * (tau2_std - tau1_std) + alpha1;
  
  // Priors
  u ~ normal(0, sigma_u);
  v ~ normal(0, sigma_v);
  sigma_u ~ exponential(2);
  sigma_v ~ exponential(2);
  alpha1 ~ normal(alpha1_mu, 2);
  beta_k ~ normal(0, 1);
  // beta_l ~ normal(0, 2);
  tau1 ~ uniform(5, 9);
  tau2 ~ uniform(10, 17);
  
  vector[N] abk; // Effect of age and breakpoints on breeding linear predictor
  vector[N] eta; // Linear predictor
  for (i in 1:N) {
    // Effect of age and breakpoints
    if (k[i] < tau1_std) {
      abk[i] = alpha1 + beta_k[1] * (k_std[i] - tau1_std);
    } else if (k[i] < tau2_std) {
      abk[i] = alpha1 + beta_k[2] * (k_std[i] - tau1_std);
    } else {
      abk[i] = alpha2 + beta_k[3] * (k_std[i] - tau2_std);
    }
    
    // Linear predictor
    // eta[i] = abk[i] + beta_l * l[i] + u[a[i]] + v[y[i]];
    eta[i] = abk[i] + u[a[i]] + v[y[i]];
    
    b[i] ~ bernoulli(inv_logit(eta[i]));
  }
}
