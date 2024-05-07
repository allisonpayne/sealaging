data {
  int<lower=0> n_ind;
  int<lower=0> n_occ;
  array[n_ind, n_occ] int<lower=0, upper=1> obs;
  array[n_ind, n_occ] int age;
}
parameters {
  // Phi
  // Fixed effects of age on survival and breeding in transition matrix
  matrix[2, 2] alpha_surv_age; // Fixed intercept of age on survival
  matrix[2, 2] beta_surv_age;  // Fixed intercept and slope of age on survival
  matrix[2, 2] alpha_bred_age; // Fixed intercept of age on breeding (conditional on survival)
  matrix[2, 2] beta_bred_age;  // Fixed slope of age on breeding (conditional on survival)
  
  // Psi
  // Random effect of year on detectability by state.
  vector[2] psi_eta[n_occ];    // detectability per year by state
  real<lower=0, upper=10> psi_sigma[2]; // std dev of psi
}
transformed parameters {
  // States:
  // 1 = breeding
  // 2 = non-breeding
  // 3 = dead
  array[n_ind, n_occ, 3] simplex[3] phi; // Transition matrix
  
  // Set transition matrix (phi) with survival and breeding
  for (i in 1:n_ind) {
    for (t in 1:n_occ) {
      // If not sexually mature yet (age <3), doesn't matter
      if (age[i, j] >= 3) {
        // Otherwise, transition probabilities have a fixed effect from age
        real surv_eta = alpha_surv_age + beta_surv_age * age[i, j];
        real bred_eta = alpha_bred_age + beta_bred_age * age[i, j];
        for (k in 1:2) {
          phi[i, t, k] = softmax([surv_eta[k, 1] * bred_eta[k, 1],       // Survived + bred
                                  bred_eta[k, 2] * (1 - bred_eta[k, 2]), // Suvived no breed
                                  1 - surv_eta[1]]);                     // Died
        }
        phi[i, t, 3] = [0, 0, 1]';  // Once dead always dead
      }
    }
  }
  
  psi = inv_logit(psi_eta)
}
model {
  array[3] real acc;
  array[n_occ] vector[3] gamma;
  
  // Psi prior
  psi_eta ~ normal(0, psi_sigma);
  
  // Phi priors
  alpha_surv_age ~ normal(0, 10);
  beta_surv_age ~ normal(0, 10);
  alpha_bred_age ~ normal(0, 10);
  beta_bred_age ~ normal(0, 10);

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1:n_ind) {
    for (t in 1:n_occ) {
      if (age[i, t] < 3) {
        for (k in 1:3) {
          gamma[t, k] = 1;
        }
      } else {
        for (j in 1:3) {
          acc[j] = gamma[t - 1, j] * ps[j, i, t - 1, k]
                     * po[k, i, t - 1, y[i, t]];
        }
      }
      for (k in 1:3) {
        
      }
      if (age[i, t] < 3) {
        for (k in 1:3) {
          gamma[t, k] = 1;
        }
      } else {
        gamma[t, i] = gamma[t - 1, ] psi[] obs[i, t]
      }
    }
    if (first[i] > 0) {
      for (k in 1 : 3) {
        gamma[first[i], k] = y[i, first[i]] == k;
      }
      
      for (t in (first[i] + 1) : n_occasions) {
        for (k in 1 : 3) {
          for (j in 1 : 3) {
            acc[j] = gamma[t - 1, j] * ps[j, i, t - 1, k]
                     * po[k, i, t - 1, y[i, t]];
          }
          gamma[t, k] = sum(acc);
        }
      }
      target += log(sum(gamma[n_occasions]));
    }
  }
}
