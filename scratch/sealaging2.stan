// States:
// 1 = breeding
// 2 = non-breeding
// 3 = dead
// Observations:
// 1 = obs breeding
// 2 = obs non-breeding
// 3 = not obs
data {
  int<lower=0> n_ind;           // Number of seals
  int<lower=0> n_occ;           // Number of occassions (years) in study as a whole (e.g., 1987-2022)
  array[n_ind, n_occ] int obs;  // Observation matrix (all 1s, 2s, and 3s)
  array[n_ind, n_occ] int age;  // Age matrix i.e., age of seal across all years. Can be negative!
}
parameters {
  // Initial state = probability of breeding at maturity (age 4)
  real<lower=0,upper=1> pb_maturity;
  
  // Logit probability of survival, breeding, and detection.
  vector[2] lps;
  vector[2] lpb;
  matrix[2, n_occ] rlpd; // random effect of year on p(detection)
  vector[2] rlpd_mu;  // mean of year random effect on p(detection)
  vector<lower=0>[2] rlpd_sigma; // dispersion of year random effect on p(detection)
}
transformed parameters {
  // invert logit tranformation to get p(survival), p(breeding), p(detection)
  vector[2] ps = inv_logit(lps);
  vector[2] pb = inv_logit(lpb);
  matrix[2, n_occ] pd = inv_logit(rlpd);
}
model {
  // We use acc and gamma to calculate the likelihood of the observations given 
  // our parameters. It's called the "forward algorithm" in HMM literature.
  array[3] real acc;
  array[n_occ] vector[3] gamma;
  // Transition and emission matrices
  array[3] vector[3] theta;  // transition matrix
  array[3, n_occ] vector[3] phi;    // emission matrix
  
  // Transitions between states are mediated by survival (ps) and transitions
  // between states 1 (breeding) and 2 (non-breeding), theta12.
  theta[1] = [ps[1] * pb[1],
              ps[1] * (1 - pb[1]),
              1 - ps[1]]';
  theta[2] = [ps[2] * pb[2],       
              ps[2] * (1 - pb[2]), 
              1 - ps[2]]';           
  theta[3] = [0, 0, 1]';
  
  // Only some emissions are possible
  // e.g., emissions for state 1 (breeding) don't allow for "observed not breeding"
  // This allows probability of detection to vary between breeding and non-breeding states
  for (t in 1:n_occ) {
    phi[1, t] = [pd[1, t], 0, 1 - pd[1, t]]';   
    phi[2, t] = [0, pd[2, t], 1 - pd[2, t]]';
    phi[3, t] = [0, 0, 1]';
  }
  
  // Priors
  // Probability of breeding at maturity (age 4) concentrated between 35% and 65%
  pb_maturity ~ beta(3, 3);
  // Probability of survival concentrated on 54-82% in both states 1 and 2
  lps ~ normal(logit(0.7), 1); // note lps is on logit scale
  // Probability of breeding concentrated on 43-75% in both states 1 and 2
  lpb ~ normal(logit(0.6), 1);
  // Random effect of year on detection
  for (t in 1:n_occ) {
    rlpd[1, t] ~ normal(rlpd_mu[1], rlpd_sigma[1]);
    rlpd[2, t] ~ normal(rlpd_mu[2], rlpd_sigma[2]);
  }
  rlpd_mu[1] ~ normal(logit(0.9), 1);
  rlpd_mu[2] ~ normal(logit(0.25), 1);
  rlpd_sigma ~ cauchy(0, 1);
  
  // When did each animal reach maturity (i.e. age 4)?
  int maturity[n_ind];
  for (i in 1:n_ind) {
    maturity[i] = 4 - age[i, 1] + 1;
  }
  
  // Likelihood
  // Calculated using forward algorithm
  for (i in 1:n_ind) {
    // Initialize at sexual maturity (age 4)
    if (obs[i, maturity[i]] == 1) {
      gamma[maturity[i]] = log([1, 0, 0])';
    } else if (obs[i, maturity[i]] == 2) {
      gamma[maturity[i]] = log([0, 1, 0])';
    } else {
      gamma[maturity[i]] = log([pb_maturity, 1 - pb_maturity, 0])';
    }
    // Progressing through life haha
    for (t in (maturity[i] + 1):n_occ) {
      for (k in 1:3) {
        for (j in 1:3) {
          // HMMs with constraints are non-trivial!
          // https://discourse.mc-stan.org/t/hidden-markov-model-with-constraints/1625/7
          if (theta[j, k] > 0 && phi[k, t, obs[i, t]] > 0) {
            // IF transition is possible (theta[j, k] > 0) AND emission is 
            // possible (phi[k, obs[i, t]] > 0) THEN update accumulator
            acc[j] = gamma[t - 1, j];
            acc[j] += log(theta[j, k]);
            acc[j] += log(phi[k, t, obs[i, t]]);
          } else {
            acc[j] = 0;
          }
        }
        gamma[t, k] = log_sum_exp(acc);
      }
    }
    target += log_sum_exp(gamma[n_occ]);;
  }
}
