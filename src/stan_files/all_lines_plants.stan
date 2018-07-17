


functions {
    int n_noninf_vector(vector x) {
        int n_noninf = 0;
        for (i in 1:rows(x)) {
            if (!is_inf(x[i])) n_noninf += 1;
        }
        return n_noninf;
    }
    vector remove_inf_vector(vector x) {
        vector[n_noninf_vector(x)] out;
        int j = 1;
        for (i in 1:rows(x)) {
            if (!is_inf(x[i])) {
                out[j] = x[i];
                j += 1;
            }
        }
        return out;
    }
    int n_nonzeros_int(int[] x) {
        int n_nonzero = 0;
        for (i in 1:(dims(x)[1])) {
            if (x[i] != 0) n_nonzero += 1;
        }
        return n_nonzero;
    }
    int[] remove_zeros_int(int[] x) {
        int out[n_nonzeros_int(x)];
        int j = 1;
        for (i in 1:(dims(x)[1])) {
            if (x[i] != 0) {
                out[j] = x[i];
                j += 1;
            }
        }
        return out;
    }
}
data {
    int<lower=1> N_ts;                              // Number of time series (line + rep)
    int<lower=1> max_nobs;                          // Max # observations per time series
    int<lower=0> max_miss;                          // Max # missing per time series
    int<lower=1, upper=max_nobs+max_miss> max_time; // Max # time points per times series
    int<lower=1> n_lines;                           // Number of aphid lines
    int<lower=0> ii_obs[N_ts, max_nobs];            // indices for observations (0 --> NA)
    int<lower=0> ii_miss[N_ts, max_miss];           // indices for missing (0 --> NA)
    int<lower=1, upper=n_lines> L[N_ts];            // aphid line for each time series
    matrix<lower=0>[max_nobs, N_ts] X_obs;          // observed log(N_t) (Inf --> NA)

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

    // Missing log(N_t):
    matrix[max_miss, N_ts] X_miss;

    // Z-transforms:

    real Z_theta;
    real Z_phi;
    real Z_s_theta;
    real Z_s_phi;
    real Z_hats_phi;
    vector[n_lines] Z_R;
    vector[n_lines] Z_A;
    vector[N_ts] Z_P;
    real Z_s_epsilon;

    matrix[max_time - 1, N_ts] Z_X;

}
transformed parameters {

    // Means (on transformed scale):
    real theta;                             // mean growth rates (r)
    real phi;                               // mean density dependences (alpha)
    // SDs (on transformed scale):
    real s_theta;                           // variation in r between lines
    real s_phi;                             // variation in alpha between lines
    real hats_phi;                          // variation in alpha within lines

    // Final estimates:
    vector[n_lines] R;                      // growth rates (r)
    vector[n_lines] A;                      // density dependences (alpha) per aphid line
    vector[N_ts] P;                         // density dependences per time series
    real s_epsilon;                         // SD of process error

    matrix[max_time, N_ts] X;      // Missing and observed log(N_t)

    theta = mu_theta + sigma_theta * Z_theta;
    phi = mu_phi + sigma_phi * Z_phi;

    s_theta = exp(gamma + sigma_gamma * Z_s_theta);
    s_phi = exp(delta + sigma_delta * Z_s_phi);
    hats_phi = exp(zeta + sigma_zeta * Z_hats_phi);

    R = exp(theta + s_theta * Z_R);
    A = inv_logit(phi + s_phi * Z_A);
    s_epsilon = exp(tau + sigma_tau * Z_s_epsilon);


    /*
     -----------
     Missing and observed log(N_t)
     -----------
    */
    for (j in 1:N_ts) {
        int n_;
        /*
         Deal with missing data
        */
        {
            // Calculate number of missing and observed time points:
            int n_obs_ = n_nonzeros_int(ii_obs[j,]);
            int n_miss_ = n_nonzeros_int(ii_miss[j,]);
            // Make vector of observed time point indices
            int obs_inds[n_obs_] = remove_zeros_int(ii_obs[j,]);
            // Make vector of observed time point log(N_t) values:
            vector[n_obs_] X_obs_ = remove_inf_vector(X_obs[,j]);
            X[obs_inds, j] = X_obs_;
            // If there are missing values, make vector of observed time point indices
            // and insert them:
            if (n_miss_ > 0) {
                int miss_inds[n_miss_] = remove_zeros_int(ii_miss[j,]);
                for (i_ in 1:n_miss_) {
                    X[miss_inds[i_], j] = X_miss[i_, j];
                }
            }
            // Needed for below:
            n_ = n_obs_ + n_miss_;
            // Fill in these cells bc I can't have any be empty:
            if (n_ < max_time) for (t in (n_+1):max_time) X[t, j] = 0;
        }
        /*
         Compute from the actual Ricker model
        */
        {
            // R for this time series:
            real r_ = R[L[j]];
            // A for this time series:
            real a_ = inv_logit(phi + s_phi * Z_A[L[j]] + hats_phi * Z_P[j]);
            P[j] = a_;
            /*
             Now filling in `X`.
             (Using for loop to avoid warning about making a copy)
            */
            for (t in 1:(n_-1)) {
                X[t+1, j] = X[t, j] + r_ * (1 - a_ * exp(X[t, j])) +
                    (Z_X[t, j] * s_epsilon);
            }
        }
    }

}
model {

    Z_theta ~ normal(0, 1);
    Z_phi ~ normal(0, 1);
    Z_s_theta ~ normal(0, 1);
    Z_s_phi ~ normal(0, 1);
    Z_hats_phi ~ normal(0, 1);
    Z_R ~ normal(0, 1);
    Z_A ~ normal(0, 1);
    Z_P ~ normal(0, 1);
    Z_s_epsilon ~ normal(0, 1);
    for (i in 1:cols(Z_X)) Z_X[,i] ~ normal(0, 1);

}
