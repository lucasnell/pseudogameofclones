functions {

    vector ricker(vector X, int n_, real r_, real a_) {
        int X_size = rows(X);
        vector[X_size] X_out;
        X_out[1] = X[1];
        X_out[2:n_] = X[1:(n_-1)] + r_ * (1 - a_ * exp(X[1:(n_-1)]));
        if (n_ < X_size) for (t in (n_+1):X_size) X_out[t] = 0;
        return X_out;
    }

}
data {

    int<lower=1> N_ts;                          // Number of time series (line + rep)
    int<lower=1> max_reps;                      // Max reps per time series
    int<lower=1> n_lines;                       // number of aphid lines
    int<lower=1> nobs_ts[N_ts];                 // # observations for each time series
    int<lower=1, upper=n_lines> L[N_ts];        // aphid line for each time series
    matrix<lower=0>[max_reps, N_ts] X;          // log(N_t)
    // Vector of priors:
    real theta[12];
    // Priors are as follows:
    //     1.  Mean for process error SD
    //     2.  SD for process error SD
    //     3.  Mean for log(growth rate) mean
    //     4.  SD for log(growth rate) mean
    //     5.  Mean for log(growth rate) SD
    //     6.  SD for log(growth rate) SD
    //     7.  Mean for logit(density dependence) mean
    //     8.  SD for logit(density dependence) mean
    //     9.  Mean for logit(density dependence) among-line SD
    //     10. SD for logit(density dependence) among-line SD
    //     11. Mean for logit(density dependence) within-line SD
    //     12. SD for logit(density dependence) within-line SD

}
parameters {

    // Z-transforms:

    vector[n_lines] Z_r;
    vector[n_lines] Z_a_a;
    vector[N_ts] Z_a_w;

    real sigma_epsilon;                 // process error
    // Means and SDs (on transformed scale):
    real rho;                           // mean growth rates: log(r)
    real sigma_rho;                     // variation in r between lines
    real phi;                           // mean density dependences: logit(alpha)
    real sigma_phi_a;                   // variation in alpha among lines
    real sigma_phi_w;                   // variation in alpha within lines

}
transformed parameters {

    matrix[max_reps, N_ts] X_pred;      // predicted X not including process error

    for (j in 1:N_ts) {
        // number of observations for this time series:
        int n_ = nobs_ts[j];
        // R for this time series:
        real r_ = exp(rho + sigma_rho * Z_r[L[j]]);
        // A for this time series:
        real a_ = inv_logit(phi + sigma_phi_a * Z_a_a[L[j]] + sigma_phi_w * Z_a_w[j]);
        // Now filling in `X_pred`:
        X_pred[, j] = ricker(X[, j], n_, r_, a_);
    }

}
model {

    Z_r ~ normal(0, 1);                 // for growth rates by line
    Z_a_a ~ normal(0, 1);               // for density dependence by line
    Z_a_w ~ normal(0, 1);               // for density dependence by plant

    sigma_epsilon ~ normal(theta[1], theta[2])T[0,];
    rho  ~ normal(theta[3], theta[4]);
    sigma_rho  ~ normal(theta[5], theta[6])T[0,];
    phi  ~ normal(theta[7], theta[8]);
    sigma_phi_a  ~ normal(theta[9], theta[10])T[0,];
    sigma_phi_w  ~ normal(theta[11], theta[12])T[0,];

    for (j in 1:N_ts) {
        X[2:nobs_ts[j], j] ~ normal(X_pred[2:nobs_ts[j], j], sigma_epsilon);
    }

}
