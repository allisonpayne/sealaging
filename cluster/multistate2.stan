//Adapted from https://github.com/stan-dev/example-models/blob/master/BPA/Ch.09/ms3_multinomlogit.stan
// -------------------------------------------------
// States (S):
// 1 Breeding
// 2 Non-Breeding
// 3 Dead
// Observations (O):
// 1 Observed breeding
// 2 Observed non-breeding
// 3 Not seen
// -------------------------------------------------
//Mark Recap with time varying detection probability 

functions {
  /**
   * Return an integer value denoting occasion of first capture.
   * This function is derived from Stan Modeling Language
   * User's Guide and Reference Manual.
   *
   * @param y         Observed values
   * @return Occasion of first capture
   */
  int first_capture(array[] int y_i) {
    for (k in 1 : size(y_i)) {
      if (y_i[k] != 3) {
        return k;
      }
    }
    return 0;
  }
  
  /**
   * Return a simplex such as follows (thanks to Bob Carpenter):
   * p[1] <- exp(lp[1]) / (1.0 + exp(lp[1]) + exp(lp[2]));
   * p[2] <- exp(lp[2]) / (1.0 + exp(lp[1]) + exp(lp[2]));
   * p[3] <- 1.0 - p[1] - p[2];
   *
   * @param lp   N-dimension vector
   * @return (N+1)-simplex of given vector and 0
   */
  vector softmax_0(vector lp) {
    vector[num_elements(lp) + 1] lp_temp;
    
    lp_temp[1 : num_elements(lp)] = lp;
    lp_temp[num_elements(lp) + 1] = 0;
    return softmax(lp_temp);
  }
}
data {
  int<lower=0> nind;
  int<lower=0> n_occasions;
  array[nind, n_occasions] int<lower=1, upper=3> y;
}
transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  array[nind] int<lower=0, upper=n_occasions> first;
  
  for (i in 1 : nind) {
    first[i] = first_capture(y[i]);
  }
}
parameters {
  real<lower=0, upper=1> phiA; // Survival probability if breeding
  real<lower=0, upper=1> phiB; // Survival probability if non-breeding
  vector<lower=0, upper=1>[n_occasions] pA; // Detection probability if breeding
  vector<lower=0, upper=1>[n_occasions] pB; // Detection probability if non-breeding
  vector[1] lpsiA; // Logit of movement probability from breeding
  vector[1] lpsiB; // Logit of movement probability from non-breeding
}
transformed parameters {
  simplex[2] psiA; // Movement probability from breeding
  simplex[2] psiB; // Movement probability from non breeding
  array[3, nind, n_occ_minus_1] simplex[3] ps;
  array[3, nind, n_occ_minus_1] simplex[3] po;
  
  // Constrain the transitions such that their sum is < 1
  psiA = softmax_0(lpsiA);
  psiB = softmax_0(lpsiB);

  // Define state-transition and observation matrices
  for (i in 1 : nind) {
    // Define probabilities of state S(t+1) given S(t)
    for (t in 1 : n_occ_minus_1) {
      ps[1, i, t, 1] = phiA * psiA[1];
      ps[1, i, t, 2] = phiA * psiA[2];
      ps[1, i, t, 3] = 1.0 - phiA;
      ps[2, i, t, 1] = phiB * psiB[1];
      ps[2, i, t, 2] = phiB * psiB[2];
      ps[2, i, t, 3] = 1.0 - phiB;
      ps[3, i, t, 1] = 0.0;
      ps[3, i, t, 2] = 0.0;
      ps[3, i, t, 3] = 1.0;

      // Define probabilities of O(t) given S(t)
      po[1, i, t, 1] = pA[t];
      po[1, i, t, 2] = 0.0;
      po[1, i, t, 3] = 1.0 - pA[t];
      po[2, i, t, 1] = 0.0;
      po[2, i, t, 2] = pB[t];
      po[2, i, t, 3] = 1.0 - pB[t];
      po[3, i, t, 1] = 0.0;
      po[3, i, t, 2] = 0.0;
      po[3, i, t, 3] = 1.0;
    }
  }
}
model {
  array[3] real acc;
  array[n_occasions] vector[3] gamma;
  
  // Priors
  // Survival and recapture: uniform
  // Uniform priors are implicitly defined.
  //  phiA ~ uniform(0, 1);
  //  phiB ~ uniform(0, 1);
  //  pA ~ uniform(0, 1);
  //  pB ~ uniform(0, 1);

  // Normal priors on logit of all but one transition probs
  lpsiA ~ normal(0, sqrt(1000));
  lpsiB ~ normal(0, sqrt(1000));

  // Likelihood
  // Forward algorithm derived from Stan Modeling Language
  // User's Guide and Reference Manual
  for (i in 1 : nind) {
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
