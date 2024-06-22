data {
  int N;          // number of obs (seal x year)
  int A;          // number of distinct animals
  int Y;          // number of distinct years

  array[N] int b; // breeding status (0 = non-breeding, 1 = breeding)
  vector[N] k;    // ages
  array[N] int a; // animal IDs
  array[N] int y; // year IDs
  vector[N] l;    // longevities
}
parameters {
  vector[A] u;           // Ranef of animal
  vector[Y] v;           // Ranef of year
  real<lower=0> sigma_u; // SD of animal ranef
  real<lower=0> sigma_v; // SD of year ranef
  real alpha;            // Intercept (i.e. p(breeding) at breakpoint)
  vector[2] beta_k;      // Age coef pre- and post-senescence
  real beta_l;           // Coefficient for longevity
  real tau;              // Breakpoint (age at senescence)
}
model {
  // Prior on alpha is the average breeding probability
  real alpha_mu = logit(mean(b));
  
  // Priors
  u ~ normal(0, sigma_u);
  v ~ normal(0, sigma_v);
  sigma_u ~ exponential(2);
  sigma_v ~ exponential(2);
  alpha ~ normal(alpha_mu, 2);
  beta_k ~ normal(0, 2);
  beta_l ~ normal(0, 2);
  tau ~ normal(0, 2);
  
  vector[N] abk; // Effect of age and breakpoint on breeding linear predictor
  vector[N] eta; // Linear predictor
  for (i in 1:N) {
    // Effect of age and breakpoint
    if (k[i] < tau) {
      abk[i] = alpha + beta_k[1] * (k[i] - tau);
    } else {
      abk[i] = alpha + beta_k[2] * (k[i] - tau);
    }
    
    // Linear predictor
    eta[i] = abk[i] + beta_l * l[i] + u[a[i]] + v[y[i]];
    
    b[i] ~ bernoulli(inv_logit(eta[i]));
  }
}
