
functions {

    vector ricker(vector X, int start, int end, real r_, real a_) {
        vector[end - start + 1] X_out;
        X_out[1] = X[start];
        X_out[2:(end - start + 1)] = X[start:(end - 1)] +
            r_ * (1 - a_ * exp(X[start:(end - 1)]));
        return X_out;
    }

    // For using nondimensionalized versions:
    vector ricker_hat(vector X, int start, int end, real r_, real a_,
                      real mu_, real tau_) {
        vector[end - start + 1] X_out;
        X_out[1] = X[start];
        X_out[2:(end - start + 1)] = X[start:(end - 1)] +
            r_ * (1 - a_ * exp(tau_ * X[start:(end - 1)] + mu_));
        return X_out;
    }
}
