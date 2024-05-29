// States:
// 1 = breeding
// 2 = non-breeding
// 3 = dead
// Observations:
// 1 = obs breeding
// 2 = obs non-breeding
// 3 = not obs
functions {
  
}
data {
  int<lower=0> n_ind;           // Number of seals
  int<lower=0> n_occ;           // Number of occassions (years) in study as a whole
  array[n_ind, n_occ] int obs;  // Observation matrix (all 1s, 2s, and 3s)
  array[n_ind, n_occ] int age;  // Age matrix i.e., age of seal across all years. Can be negative!
}
parameters {
  // Initial state = probability of breeding at maturity (age 4)
  real<lower=0,upper=1> pb_maturity;
  
  // Logit-scale probability of survival, breeding, and detection.
  vector[2] logit_ps;
  vector[2] logit_pb;
  matrix[2, n_occ] logit_pd;         // random effect of year on p(detection)
  vector[2] logit_pd_mu;             // mean of year random effect on p(detection)
  vector<lower=0>[2] logit_pd_sigma; // sd of year random effect on p(detection)
}
model {
  // We use acc and gamma to calculate the likelihood of the observations given 
  // our parameters. It's called the "forward algorithm" in HMM literature.
  vector[3] acc;
  array[n_occ] vector[3] gamma;
  
  // Log probabilities
  vector[2] lps = log(inv_logit(logit_ps));        // Survival
  vector[2] lpm = log1m_exp(lps);                  // Mortality
  vector[2] lpb = log(inv_logit(logit_pb));        // Breeding
  matrix[2, 2] ltheta12 = [lpb', log1m_exp(lpb)']; // Moving between states 1, 2
  matrix[2, n_occ] lpd = log(inv_logit(logit_pd)); // Detection
  
  // Priors
  // Probability of breeding at maturity (age 4) gets a flat prior
  pb_maturity ~ beta(1, 1);
  // Probability of survival concentrated on 54-82% in both states 1 and 2
  logit_ps ~ normal(logit(0.7), 1); // note lps is on logit scale
  // Probability of breeding concentrated on 43-75% in both states 1 and 2
  logit_pb ~ normal(logit(0.6), 1);
  // Random effect of year on detection
  logit_pd_mu[1] ~ normal(logit(0.9), 1);
  logit_pd_mu[2] ~ normal(logit(0.4), 1);
  logit_pd_sigma ~ cauchy(0, 5);
  logit_pd[1] ~ normal(logit_pd_mu[1], logit_pd_sigma[1]);
  logit_pd[2] ~ normal(logit_pd_mu[2], logit_pd_sigma[2]);
             
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
      gamma[maturity[i]] = [0, negative_infinity(), negative_infinity()]';
    } else if (obs[i, maturity[i]] == 2) {
      gamma[maturity[i]] = [negative_infinity(), 0, negative_infinity()]';
    } else {
      gamma[maturity[i]] = log([pb_maturity, 1 - pb_maturity, 0])';
    }
    // Likelihood of observations given parameters
    for (t in (maturity[i] + 1):n_occ) {
      if (obs[i, t] < 3) {
        // If state known (obs = 1 or 2)...
        for (j in 1:2) {
          acc[j] = gamma[t - 1, j];         // Start with previous time step
          acc[j] += lps[j];                 // Did the animal survive?
          acc[j] += ltheta12[j, obs[i, t]]; // If so, did it breed?
          acc[j] += lpd[obs[i, t], t];      // If so, did we observe it?
        }
        acc[3] = negative_infinity(); 
        gamma[t, obs[i, t]] = log_sum_exp(acc);        // Update gamma for known state
        gamma[t, 3 - obs[i, t]] = negative_infinity(); // Gamma for other states
        gamma[t, 3] = negative_infinity();             // Can't be dead
      } else {
        // If state unknown (obs = 3)
        // For transitions *to* states 1 and 2...
        for (k in 1:2) {
          for (j in 1:2) {
            acc[j] = gamma[t - 1, j];     // Start with previous time step
            acc[j] += lps[j];             // Did the animal survive?
            acc[j] += ltheta12[j, k];     // If so, did it breed?
            acc[j] += lpd[k, t];          // If so, did we observe it?
          }
          acc[3] = negative_infinity();   // can't transition from dead
          gamma[t, k] = log_sum_exp(acc);
        }
        // For transition *to* state 3
        for (j in 1:2) {
          acc[j] = gamma[t - 1, j] + lpm[j];
        }
        acc[3] = gamma[t - 1, 3];
        gamma[t, 3] = log_sum_exp(acc);
      }
    }
    target += log_sum_exp(gamma[n_occ]);;
  }
}
