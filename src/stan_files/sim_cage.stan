
functions {

    int find_int_array(int x, int[] y) {
        int result = 0;
        for (i in 1:(dims(y)[1])) {
            if (y[i] == x) result = 1;
        }
        return result;
    }

    int get_n_replaced(int[,] repl_plants, int repl_ind) {

        int n_repl_plants = 0;
        for (i_ in 1:(dims(repl_plants)[2])) {
            int i = repl_plants[repl_ind, i_];
            if (i != 0) n_repl_plants += 1;
        }
        return n_repl_plants;
    }

    int[] get_not_replaced(int n_non_repl_plants, int n_plants,
                           int[,] repl_plants, int repl_ind) {
        int non_repl_plants[n_non_repl_plants] = rep_array(0, n_non_repl_plants);
        int ii = 1;
        for (i in 1:n_plants) {
            if (!find_int_array(i, repl_plants[repl_ind,])) {
                non_repl_plants[ii] = i;
                ii += 1;
            }
        }
        return non_repl_plants;
    }
}
data {
    int<lower=1> n_plants;                  // Number of plants per cage
    int<lower=1> n_lines;                   // number of aphid lines per cage
    matrix<lower=0>[n_plants, n_lines] X_0; // log(N) at time t=0
    int<lower=1> max_t;                     // Time steps per cage
    vector<lower=0>[n_lines] R;             // Max growth rates per aphid line
    vector<lower=0,upper=1>[n_lines] A;     // Density dependence per aphid line
    vector[n_lines] D_slope;                // Dispersal slope per aphid line
    vector[n_lines] D_inter;                // Dispersal intercept per aphid line
    real<lower=0> process_error;            // SD of process error
    vector[2] plant_mort_coefs;             // Coefficients (b0, then b1) for aphid
                                            // mortality after plant starts dying
    int<lower=1> plant_death_age;           // # days after which a plant starts dying
    int<lower=1> n_repl;                    // Number of plant replacements
    int<lower=1, upper=n_plants> max_repl;  // Max # plants replaced per replacement
    int<lower=0> repl_times[n_repl];        // Time points when plant replacements occur
    // Plants to replace when replacements occur (zeros are ignored):
    int<lower=0, upper=n_plants> repl_plants[n_repl, max_repl];
}
generated quantities {

    matrix[max_t+1, n_lines] X_out[n_plants];     // 3D array for output

    // Fill with initial values and -Inf for all other rows:
    for (i in 1:n_plants) {
        X_out[i] = rep_matrix(negative_infinity(), max_t+1, n_lines);
        for (j in 1:n_lines) {
            X_out[i][1,j] = X_0[i,j];
        }
    }

    {
        // Matrix to keep track of extinctions
        matrix[n_plants, n_lines] extinct = rep_matrix(0, n_plants, n_lines);
        // Vector to keep track of plant ages:
        int plant_ages[n_plants] = rep_array(0, n_plants);
        // Index to keep track of position in `repl_times`
        int repl_ind = 1;
        // Matrices keeping track of numbers of dispersed aphids:
        matrix[n_lines, n_plants] emigrants = rep_matrix(0, n_lines, n_plants);
        matrix[n_lines, n_plants] immigrants = rep_matrix(0, n_lines, n_plants);
        // mortality at a time x after plant death starts, to keep from having to do
        // this calculation many times. `max_t` should be more than enough.
        vector[max_t] morts;
        for (x in 1:max_t) {
            real x_ = x;
            morts[x] = inv_logit(plant_mort_coefs[1] + plant_mort_coefs[2] * x_);
        }

        for (t in 1:max_t) {

            /*
             Go through once to get/update some parameters first:
            */
            // Dispersal lambdas (for Poisson distr.) for all lines and plants
            matrix[n_lines, n_plants] D_lambdas;
            // Summed density dependences * N_t among all lines for each plant
            vector[n_plants] Z = rep_vector(0, n_plants);
            for (i in 1:n_plants) {
                // Update plant ages to time t+1:
                plant_ages[i] += 1;
                for (j in 1:n_lines) {
                    D_lambdas[j,i] = exp(D_inter[j] + D_slope[j] * X_out[i][t,j]);
                    Z[i] += (A[j] * exp(X_out[i][t,j]));
                }
            }
            // Generate numbers of dispersed aphids:
            for (j in 1:n_lines) {
                int dispersed_;
                for (from_i in 1:n_plants) {
                    for (to_i in 1:n_plants) {
                        if (from_i == to_i) continue;
                        dispersed_ = poisson_rng(D_lambdas[j,from_i] / (n_plants - 1));
                        emigrants[j, from_i] += dispersed_;
                        immigrants[j, to_i] += dispersed_;
                    }
                }
            }

            /*
             Go back through to simulate at time t+1
            */
            for (i in 1:n_plants) {

                for (j in 1:n_lines) {

                    // Calculate the net influx of this aphid line from other plants.
                    real immigration = immigrants[j,i];
                    real emigration = emigrants[j,i];

                    // If it's extinct and no one's coming in, skip the rest
                    // if (extinct[i,j] == 1 && immigration == 0) continue;

                    // Start out new X based on previous time step
                    X_out[i][t+1,j] = X_out[i][t,j];
                    /*
                     Now add growth and density dependence.
                     If the plant is past plant-death age, then plant-death-induced
                     mortality overrides replaces normal growth and density dependence,
                     and will be added after adding dispersal.
                    */
                    if (plant_ages[i] <= plant_death_age) {
                        X_out[i][t+1,j] += (R[j] * (1 - Z[i]));
                    }
                    // Add process error:
                    X_out[i][t+1,j] += (normal_rng(0, 1) * process_error);

                    // Temporarily convert X to units of N, to add dispersal:
                    X_out[i][t+1,j] = exp(X_out[i][t+1,j]) + immigration - emigration;

                    /*
                     If the plant is past plant-death age, then add
                     plant-death-induced mortality here.
                    */
                    if (plant_ages[i] > plant_death_age) {
                        X_out[i][t+1,j] *= morts[(plant_ages[i] - plant_death_age)];
                    }

                    // Convert back to units of X, but check for extinction first
                    if (X_out[i][t+1,j] < 1) {
                        X_out[i][t+1,j] = negative_infinity();
                        extinct[i,j] = 1;
                    } else X_out[i][t+1,j] = log(X_out[i][t+1,j]);

                }

            }

            /*
             If it's a replacement time point, disperse all aphids from replaced
             plants to others
            */
            if (t == repl_times[repl_ind]) {
                // Get # of plants that are being replaced (the # non-zero entries)
                int n_repl_plants = get_n_replaced(repl_plants, repl_ind);
                int n_non_repl_plants = n_plants - n_repl_plants;
                real n_non_repl_plants_ = n_non_repl_plants;
                // Fill non-replaced plant indices:
                int non_repl_plants[n_non_repl_plants] = get_not_replaced(
                    n_non_repl_plants, n_plants, repl_plants, repl_ind);

                // Update replaced plants' ages:
                for (i_ in 1:(dims(repl_plants)[2])) {
                    int i = repl_plants[repl_ind, i_];
                    if (i == 0) continue;
                    plant_ages[i] = 0;
                }

                // Adjust pools of aphid lines in replaced and non-replaced plants:
                for (j in 1:n_lines) {
                    // Pool of aphids of line j to disperse across non-replaced plants:
                    real aphid_pool_j = 0;
                    // Calculate aphid_pool_j and remove from replaced plant(s):
                    for (i_ in 1:max_repl) {
                        int i = repl_plants[repl_ind, i_];
                        if (i == 0) continue;
                        aphid_pool_j += exp(X_out[i][t+1,j]);
                        X_out[i][t+1,j] = negative_infinity();
                        extinct[i,j] = 1;
                    }
                    // Now add a portion to non-replaced plant(s):
                    aphid_pool_j /= n_non_repl_plants_;
                    for (i_ in 1:n_non_repl_plants) {
                        int i = non_repl_plants[i_];
                        // Check that it's not adding < 1 aphid to an extinct plant:
                        if (extinct[i,j] == 1 && aphid_pool_j < 1) continue;
                        // If not, add `aphid_pool_j`:
                        X_out[i][t+1,j] = log(exp(X_out[i][t+1,j]) + aphid_pool_j);
                    }
                }
                // Now adjust index for `repl_times`
                repl_ind += 1;
            }


        }

    }

}


