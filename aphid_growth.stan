data {
    int<lower=0> N;             // number of observations
    vector[N] X;                  // log(N_(t+1))
    vector[N] X0;                 // log(N_(t))
    // int<lower=0> N_ts;          // number of time series (aphid line + rep)
    // int<lower=0> reps_ts[J];    // reps per time series
    // int<lower=0> line_ts[J];    // aphid line for each time series
    // int<lower=0> n_lines;       // number of aphid lines
}
parameters {
    real R0; // growth rate
    real alpha; // density dependence
    
  // real R0[n_lines];             // estimated growth rates
  // real alpha[n_lines]; // estimated density dependences
}
transformed parameters {
  vector[N_ts] X_pred;
  X_pred = X0 + R0 * (1 - alpha * exp(X0));
}
model {
  X ~ normal(X_pred, sigma);
  
  
  
}
