//
// data {
//     int<lower=1> N_ts;                          // number of time series (line + rep)
//     int<lower=1> n_lines;                       // number of aphid lines
//     int<lower=1> reps_ts[N_ts];                 // reps per time series
//     int<lower=1, upper=n_lines> line_ts[N_ts];  // aphid line for each time series
//     ragged_vector[reps_ts] X;                   // log(N_t)
// }
// parameters {
//     real R0[n_lines];               // estimated growth rates
//     real alpha[n_lines];            // estimated density dependences
//     real sigma_process;             // process error
// }
// transformed parameters {
//     R0 = R0_mean + sigma_R0 * Z_R0;
//     alpha = alpha_mean + sigma_alpha * Z_alpha;
//     ragged_vector[reps_ts] X_pred;
//     for (i in 1:N_ts) {
//         X_pred[i][2:rep_ts[i]] = X[i][1:(rep_ts[i]-1)] + R0[line_ts[i]] *
//             (1 - alpha[line_ts[i]] * exp(X[i][1:(rep_ts[i]-1)]));
//     }
// }
// model {
//
//     // Inputs here are your priors:
//     real R0_mean ~ normal(0.29, );
//     real alpha_mean ~ normal();
//     real<lower=0> sigma_R0 ~ normal(0.0227, );
//     real<lower=0> sigma_alpha ~ normal();
//
//     real Z_R0[n_lines] ~ normal(0,1);
//     real Z_alpha[n_lines] ~ normal(0,1);
//
//     for (i in 1:N_ts) {
//         X[i][2:rep_ts[i]] ~ normal(X_pred[i][2:rep_ts[i]], sigma_process);
//     }
//
// }




data {
    int<lower=1> N;     // Number of observations
    vector[N] X;        // log(N_t)
}
parameters {
    real<lower=0> R0;               // estimated growth rate
    real<lower=0> alpha;            // estimated density dependence
    real<lower=0> sigma_process;    // process error
    // real<lower=0> R0_mean;
    // real<lower=0> alpha_mean;
    // real<lower=0> sigma_process_mean;
    // real<lower=0> sigma_R0;
    // real<lower=0> sigma_alpha;
    // real<lower=0> sigma_sigma_process;
}
transformed parameters {
    // real Z_R0;
    // real Z_alpha;
    // real Z_sigma_process;
    // R0 = R0_mean + sigma_R0 * Z_R0;
    // alpha = alpha_mean + sigma_alpha * Z_alpha;
    // sigma_process = sigma_process_mean + sigma_sigma_process * Z_sigma_process;
    vector[N] X_pred;
    X_pred[1] = X[1];
    X_pred[2:N] = X[1:(N-1)] + R0 * (1 - alpha * exp(X[1:(N-1)]));
}
model {

    // // Inputs here are your priors:
    // R0_mean ~ normal(0.29, 10);
    // alpha_mean ~ normal(509, 500);
    // sigma_process_mean ~ normal(2, 10);
    // sigma_R0 ~ normal(0.0227, 1);
    // sigma_alpha ~ normal(100, 50);
    // sigma_sigma_process ~ normal(2, 10);

    // Z_R0 ~ normal(0,1);
    // Z_alpha ~ normal(0,1);
    // Z_sigma_process ~ normal(0,1);

    R0 ~ normal(0.29, 0.0227*10);
    alpha ~ normal(509, 100 * 10);
    sigma_process ~ normal(2, 10);

    X[2:N] ~ normal(X_pred[2:N], sigma_process);

}

