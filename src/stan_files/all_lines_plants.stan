
/*
    This version differs from `all_lines.stan` because it also models a parameter
    for each plant's affect on the line's alpha parameter.
*/

data {
    int<lower=1> N_ts;                          // Number of time series (line + rep)
    int<lower=1> max_reps;                      // Max reps per time series
    int<lower=1> n_lines;                       // number of aphid lines
    int<lower=1> nobs_ts[N_ts];                 // # observations for each time series
    int<lower=1, upper=n_lines> line_ts[N_ts];  // aphid line for each time series
    matrix[max_reps, N_ts] X;                   // log(N_t)
}
parameters {
    // Means:
    real<lower=0> mean_R0;
    real<lower=0> mean_alpha;
    // real mean_beta;              // I'm going to force this to be one
    real<lower=0> mean_process;
    // SDs:
    real<lower=0> sigma_R0;
    real<lower=0> sigma_alpha;
    real<lower=0> sigma_beta;
    real<lower=0> sigma_process;
    // Z-transforms:
    vector<lower=0>[n_lines] Z_R0;
    vector<lower=0>[n_lines] Z_alpha;
    vector[N_ts] Z_beta;
    real<lower=0> Z_process;
}
transformed parameters {

    vector<lower=0>[n_lines] R0;    // estimated growth rate
    vector<lower=0>[n_lines] alpha; // estimated density dependence
    vector<lower=0>[N_ts] beta;     // effects of plants on alpha
    real<lower=0> process;          // process error
    matrix[max_reps, N_ts] X_pred;

    R0 = mean_R0 + sigma_R0 * Z_R0;
    alpha = mean_alpha + sigma_alpha * Z_alpha;
    beta = 1 + sigma_beta * Z_beta;
    process = mean_process + sigma_process * Z_process;

    for (i in 1:N_ts) {
        X_pred[1, i] = X[1, i];
        X_pred[2:nobs_ts[i], i] = X[1:(nobs_ts[i]-1), i] + R0[line_ts[i]] *
            (1 - beta[i] * alpha[line_ts[i]] * exp(X[1:(nobs_ts[i]-1), i]));
        if (nobs_ts[i] < max_reps) {
            for (j in (nobs_ts[i]+1):max_reps) X_pred[j, i] = 0;
        }
    }

}
model {

    // Inputs here are your priors:
    // sigma_process ~ normal(0.2, 3);
    // Means
    mean_R0 ~ normal(0.29, 1);
    mean_alpha ~ normal(0.0029, 1);
    // mean_beta ~ normal(0, 1);
    mean_process ~ normal(0.1, 0.000364 * 10);
    // SDs
    sigma_R0 ~ cauchy(0, 0.005);
    sigma_alpha ~ cauchy(0, 0.0005);
    sigma_beta ~ cauchy(0, 0.0005);
    sigma_process ~ cauchy(0, 0.05);
    // Z-transforms
    Z_R0 ~ normal(0,1);
    Z_alpha ~ normal(0,1);
    Z_beta ~ normal(0,1);
    Z_process ~ normal(0,1);

    for (i in 1:N_ts) {
        X[2:nobs_ts[i], i] ~ normal(X_pred[2:nobs_ts[i], i], process);
    }

}
