
data {
    int<lower=1> N;     // Number of observations
    vector[N] X;        // log(N_t)
}
parameters {
    real<lower=0> process;    // process error
    // Means:
    real<lower=0> mean_R0;
    real<lower=0> mean_alpha;
    // real<lower=0> mean_process;
    // SDs:
    real<lower=0> sigma_R0;
    real<lower=0> sigma_alpha;
    // real<lower=0> sigma_process;
    // Z-transforms:
    real<lower=0> Z_R0;
    real<lower=0> Z_alpha;
    // real<lower=0> Z_process;
}
transformed parameters {

    real<lower=0> R0;               // estimated growth rate
    real<lower=0> alpha;            // estimated density dependence
    // real<lower=0> process;          // process error
    vector[N] X_pred;

    R0 = mean_R0 + sigma_R0 * Z_R0;
    alpha = mean_alpha + sigma_alpha * Z_alpha;
    // process = mean_process + sigma_process * Z_process;

    X_pred[1] = X[1];
    X_pred[2:N] = X[1:(N-1)] + R0 * (1 - alpha * exp(X[1:(N-1)]));
}
model {

    // Inputs here are your priors:
    // sigma_process ~ normal(0.2, 3);
    // Means
    mean_R0 ~ normal(0.29, 0.0226 * 10);
    mean_alpha ~ normal(0.0029, 0.000364 * 10);
    // SDs
    sigma_R0 ~ normal(0, 1);
    sigma_alpha ~ normal(0, 0.000364 * 10);
    // Z-transforms
    Z_R0 ~ normal(0,1);
    Z_alpha ~ normal(0,1);

    // R0 ~ normal(0.29, 0.0227 * 5);
    // alpha ~ normal(0, 1);
    process ~ normal(1, 5);

    X[2:N] ~ normal(X_pred[2:N], process);

}



/*
--------
From Stan documentation for how to use Z-scores:
(Top is without z-scores, bottom is with)
--------
*/

// data {
//   int<lower=0> J;
//   real y[J];
//   real<lower=0> sigma[J];
// }
//
// parameters {
//   real mu;
//   real<lower=0> tau;
//   real theta[J];
// }
//
// model {
//   mu ~ normal(0, 5);
//   tau ~ cauchy(0, 5);
//   theta ~ normal(mu, tau);
//   y ~ normal(theta, sigma);
// }
//
//
// data {
//   int<lower=0> J;
//   real y[J];
//   real<lower=0> sigma[J];
// }
//
// parameters {
//   real mu;
//   real<lower=0> tau;
//   real theta_tilde[J];
// }
//
// transformed parameters {
//   real theta[J];
//   for (j in 1:J)
//     theta[j] = mu + tau * theta_tilde[j];
// }
//
// model {
//   mu ~ normal(0, 5);
//   tau ~ cauchy(0, 5);
//   theta_tilde ~ normal(0, 1);
//   y ~ normal(theta, sigma);
// }

