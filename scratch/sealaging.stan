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
  
  // Probability of survival, breeding, and detection.
  vector<lower=0,upper=1>[2] ps;
  vector<lower=0,upper=1>[2] pb;
  vector<lower=0,upper=1>[2] pd;
}
model {
  // We use acc and gamma to calculate the likelihood of the observations given 
  // our parameters. It's called the "forward algorithm" in HMM literature.
  array[3] real acc;
  array[n_occ] vector[3] gamma;
  // Transition and emission matrices
  array[3] vector[3] theta;  // transition matrix
  array[3] vector[3] phi;    // emission matrix
  
  // Transitions between states are mediated by survival (ps) and transitions
  // between states 1 (breeding) and 2 (non-breeding), theta12.
  theta[1] = [ps[1] * pb[1],       // breed again = p(surv|bred) * p(bred->bred)
              ps[1] * (1 - pb[1]), // skip breed = p(surv|bred) * (1-p(bred->bred)
              1 - ps[1]]';         // died = 1 - p(surv|bred)
  theta[2] = [ps[2] * pb[2],       // as above, but for state 2
              ps[2] * (1 - pb[2]), 
              1 - ps[2]]';           
  theta[3] = [0, 0, 1]';             // if you're dead you stay dead (no transitioning back)
  
  // Only some emissions are possible
  // e.g., emissions for state 1 (breeding) don't allow for "observed not breeding"
  // This allows probability of detection to vary between breeding and non-breeding states
  phi[1] = [pd[1], 0, 1 - pd[1]]';   
  phi[2] = [0, pd[2], 1 - pd[2]]';
  phi[3] = [0, 0, 1]';
  
  // Priors
  // Probability of breeding at maturity (age 4) concentrated between 35% and 65%
  pb_maturity ~ beta(3, 3);
  // Probability of survival concentrated on 65%-75% in both states 1 and 2
  ps ~ beta(6, 3);
  // Probability of breeding concentrated on 50-70% in both states 1 and 2
  pb ~ beta(6, 4);
  // Probability of detection very high in state 1, ~25% in state 2
  pd[1] ~ beta(5, 1);
  pd[2] ~ beta(2, 4);
  
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
          if (theta[j, k] > 0 && phi[k, obs[i, t]] > 0) {
            // IF transition is possible (theta[j, k] > 0) AND emission is 
            // possible (phi[k, obs[i, t]] > 0) THEN update accumulator
            acc[j] = gamma[t - 1, j] + log(theta[j, k]) 
                   + log(phi[k, obs[i, t]]);
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
