/*
 This model does NOT include among-line variance in density dependence.
 It DOES include...
 - Among-line variability in growth rates
 - Among-PLANT variability in density dependences
 - Process error that is estimated by the model
 */

#include /chunks/helpers.stan
data {

    // Indices:
    int<lower=1> n_ts;                          // number of time series (line + rep)
    int<lower=1> n_obs;                         // total number of observations
    int<lower=1> n_per[n_ts];                   // number of obs. for each time series

    // Data:
    vector<lower=0>[n_obs] X;                   // log(N_t)

    int<lower=1> n_lines;                       // number of aphid lines
    int<lower=1, upper=n_lines> L[n_ts];        // aphid line for each time series

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
    //     11. Mean for logit(density dependence) among-plant SD
    //     12. SD for logit(density dependence) among-plant SD

}
transformed data {
    vector[n_obs] X_hat;
    real mu;
    real tau;

    mu = mean(X);
    tau = sd(X);

    X_hat = (X - mu) / tau;
}
parameters {

    // Z-scores:

    vector[n_lines] Z_r;
    // vector[n_lines] Z_alpha_l;
    vector[n_ts] Z_alpha_p;

    real<lower=0> sigma_epsilon;        // process error
    // Means and SDs (on transformed scale):
    real rho;                           // mean growth rates: log(r)
    real<lower=0> sigma_rho;            // among-line SD in log(r)
    real phi;                           // mean density dependences: logit(alpha)
    // real<lower=0> sigma_phi_l;          // among-line SD in logit(alpha)
    real<lower=0> sigma_phi_p;          // within-line SD in logit(alpha)

}
transformed parameters {

    vector[n_obs] X_hat_pred;               // predicted X_hat not including process error

    // iterate over each time series
    {
        int start = 1; // starting position in `X` and `X_pred` vectors
        for (j in 1:n_ts) {
            // number of observations for this time series:
            int n_ = n_per[j];
            int end = start + n_ - 1;
            // growth rate (r) for this time series:
            real r_hat = exp(rho + sigma_rho * Z_r[L[j]]);
            // density dependence (alpha) for this time series:
            real alpha_hat = inv_logit(phi +
                // sigma_phi_l * Z_alpha_l[L[j]] +
                sigma_phi_p * Z_alpha_p[j]);
            // filling in predicted X_t+1 based on X_t:
            X_hat_pred[start:end] = ricker(X_hat, start, end, r_hat, alpha_hat);
            // iterate `start` index:
            start += n_;
        }
    }

}
model {

    Z_r ~ normal(0, 1);                 // for growth rates by line
    // Z_alpha_l ~ normal(0, 1);          // for density dependence by line
    Z_alpha_p ~ normal(0, 1);           // for density dependence by plant

    sigma_epsilon ~ normal(theta[1], theta[2])T[0,];
    rho  ~ normal(theta[3], theta[4]);
    sigma_rho  ~ normal(theta[5], theta[6])T[0,];
    phi  ~ normal(theta[7], theta[8]);
    // sigma_phi_l  ~ normal(theta[9], theta[10])T[0,];
    sigma_phi_p  ~ normal(theta[11], theta[12])T[0,];

    // Process error:
    X_hat ~ normal(X_hat_pred, sigma_epsilon);

}
generated quantities {

    // Residuals
    vector[n_obs] X_resid;

    /*
     Parameter estimates in original scale:
     */
    vector[n_obs] X_pred;
    vector[n_lines] R;  // growth rate by line
    // vector[n_lines] A;  // alpha by line
    real alpha;         // mean alpha over all lines
    real sigma_process; // SD of process error

    X_resid = X_hat - X_hat_pred;

    X_pred = X_hat_pred * tau + mu;
    R = exp(rho + sigma_rho * Z_r) * tau;
    // A = inv_logit(phi + sigma_phi_l * Z_alpha_l);
    alpha = inv_logit(phi);
    sigma_process = sigma_epsilon * tau;

}
