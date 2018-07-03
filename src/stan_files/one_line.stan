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
//         X_pred[i][2:reps_ts[i]] = X[i][1:(reps_ts[i]-1)] + R0[line_ts[i]] *
//             (1 - alpha[line_ts[i]] * exp(X[i][1:(reps_ts[i]-1)]));
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
//         X[i][2:reps_ts[i]] ~ normal(X_pred[i][2:reps_ts[i]], sigma_process);
//     }
//
// }




data {
    int<lower=1> N_ts;              // Number of time series
    int<lower=1> reps_ts[N_ts];     // reps per time series
    int<lower=1> max_reps;          // Max reps per time series
    matrix[max_reps, N_ts] X;   // log(N_t)
}
parameters {
    // real<lower=0> process;    // process error
    // Means:
    real<lower=0> mean_R0;
    real<lower=0> mean_alpha;
    real<lower=0> mean_process;
    // SDs:
    real<lower=0> sigma_R0;
    real<lower=0> sigma_alpha;
    real<lower=0> sigma_process;
    // Z-transforms:
    real<lower=0> Z_R0;
    real<lower=0> Z_alpha;
    real<lower=0> Z_process;
}
transformed parameters {

    real<lower=0> R0;               // estimated growth rate
    real<lower=0> alpha;            // estimated density dependence
    real<lower=0> process;          // process error
    matrix[max_reps, N_ts] X_pred;

    R0 = mean_R0 + sigma_R0 * Z_R0;
    alpha = mean_alpha + sigma_alpha * Z_alpha;
    process = mean_process + sigma_process * Z_process;

    for (i in 1:N_ts) {
        X_pred[1, i] = X[1, i];
        X_pred[2:reps_ts[i], i] = X[1:(reps_ts[i]-1), i] + R0 *
            (1 - alpha * exp(X[1:(reps_ts[i]-1), i]));
        if (reps_ts[i] < max_reps) {
            for (j in (reps_ts[i]+1):max_reps) X_pred[j, i] = 0;
        }
    }

}
model {

    // Inputs here are your priors:
    // sigma_process ~ normal(0.2, 3);
    // Means
    mean_R0 ~ normal(0.29, 0.0226 * 10);
    mean_alpha ~ normal(0.0029, 0.000364 * 10);
    mean_process ~ normal(0.0029, 0.000364 * 10);
    // SDs
    sigma_R0 ~ normal(0, 1);
    sigma_alpha ~ normal(0, 0.000364 * 10);
    sigma_process ~ normal(0, 0.000364 * 10);
    // Z-transforms
    Z_R0 ~ normal(0,1);
    Z_alpha ~ normal(0,1);
    Z_process ~ normal(0,1);

    // R0 ~ normal(0.29, 0.0227 * 5);
    // alpha ~ normal(0, 1);
    // process ~ normal(1, 5);

    for (i in 1:N_ts) {
        X[2:reps_ts[i], i] ~ normal(X_pred[2:reps_ts[i], i], process);
    }

}
