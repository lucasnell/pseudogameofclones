
/*
    This version differs from `all_lines.stan` because it also models variability
    in density dependence (variable `A`) within clonal lines.
*/

data {
    int<lower=1> N_ts;                          // Number of time series (line + rep)
    int<lower=1> max_reps;                      // Max reps per time series
    int<lower=1> n_lines;                       // number of aphid lines
    int<lower=1> nobs_ts[N_ts];                 // # observations for each time series
    int<lower=1, upper=n_lines> L[N_ts];        // aphid line for each time series
    matrix<lower=0>[max_reps, N_ts] X;          // log(N_t)
    // Priors for process error:
    real tau;
    real sigma_tau;
    // Priors for growth rates:
    real mu_theta;
    real sigma_theta;
    real gamma;
    real sigma_gamma;
    // Priors for density dependence:
    real mu_phi;
    real sigma_phi;
    real delta;
    real sigma_delta;
    real zeta;
    real sigma_zeta;
}
parameters {

    // Z-transforms:

    real Z_theta;
    real Z_phi;
    real Z_s_theta;
    real Z_s_phi;
    vector[n_lines] Z_S;
    vector[n_lines] Z_R;
    vector[n_lines] Z_A;
    vector[N_ts] Z_P;
    real Z_s_epsilon;

}
transformed parameters {

    // Means (on transformed scale):
    real theta;                             // mean growth rates (r)
    real phi;                               // mean density dependences (alpha)
    // SDs (on transformed scale):
    real s_theta;                           // variation in r between lines
    real s_phi;                             // variation in alpha between lines
    vector[n_lines] S;                      // variation in alpha within lines

    // Final estimates:
    vector[n_lines] R;                      // growth rates (r)
    vector[n_lines] A;                      // density dependences (alpha) per aphid line
    vector[N_ts] P;                         // density dependences per time series
    real s_epsilon;                         // SD of process error

    matrix[max_reps, N_ts] X_pred;          // predicted X not including process error


    theta = mu_theta + sigma_theta * Z_theta;
    phi = mu_phi + sigma_phi * Z_phi;

    s_theta = exp(gamma + sigma_gamma * Z_s_theta);
    s_phi = exp(delta + sigma_delta * Z_s_phi);
    S = exp(zeta + sigma_zeta * Z_S);

    R = exp(theta + s_theta * Z_R);
    A = inv_logit(phi + s_phi * Z_A);
    s_epsilon = exp(tau + sigma_tau * Z_s_epsilon);

    for (j in 1:N_ts) {
        // number of observations for this time series:
        int n_ = nobs_ts[j];
        // R for this time series:
        real r_ = R[L[j]];
        // A for this time series:
        real a_ = inv_logit(phi + s_phi * Z_A[L[j]] + S[L[j]] * Z_P[j]);
        P[j] = a_;
        // Now filling in `X_pred`:
        X_pred[1, j] = X[1, j];
        X_pred[2:n_, j] = X[1:(n_-1), j] + r_ * (1 - a_ * exp(X[1:(n_-1), j]));
        if (n_ < max_reps) for (t in (n_+1):max_reps) X_pred[t, j] = 0;
    }

}
model {

    Z_theta ~ normal(0, 1);
    Z_phi ~ normal(0, 1);
    Z_s_theta ~ normal(0, 1);
    Z_s_phi ~ normal(0, 1);
    Z_S ~ normal(0, 1);
    Z_R ~ normal(0, 1);
    Z_A ~ normal(0, 1);
    Z_P ~ normal(0, 1);
    Z_s_epsilon ~ normal(0, 1);

    for (j in 1:N_ts) {
        X[2:nobs_ts[j], j] ~ normal(X_pred[2:nobs_ts[j], j], s_epsilon);
    }

}
