
#include /chunks/helpers.stan
data {
    int<lower=1> n_plants;                  // Number of plants per cage
    int<lower=1> n_lines;                   // number of aphid lines per cage
    matrix<lower=0>[n_plants, n_lines] N_0; // N at time t=0
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

    matrix[max_t+1, n_lines] N_out[n_plants];     // 3D array for output

    // Fill with initial values and 0 for all other rows:
    for (i in 1:n_plants) {
        N_out[i] = rep_matrix(0, max_t+1, n_lines);
        for (j in 1:n_lines) {
            N_out[i][1,j] = N_0[i,j];
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
        // log(mortality) at a time x after plant death starts, to keep from having to do
        // this calculation many times. `max_t - plant_death_age` should be more
        // than enough.
        vector[(max_t - plant_death_age)] log_morts;
        for (x in 1:(max_t - plant_death_age)) {
            real x_ = x;
            log_morts[x] = log(inv_logit(plant_mort_coefs[1] + plant_mort_coefs[2] * x_));
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
                    D_lambdas[j,i] = exp(D_inter[j] + D_slope[j] * log(N_out[i][t,j]));
                    Z[i] += (A[j] * N_out[i][t,j]);
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
                    if (extinct[i,j] == 1 && immigration == 0) continue;

                    // Start out new N based on previous time step
                    N_out[i][t+1,j] = N_out[i][t,j];
                    /*
                     Now add growth, density dependence, and process error.
                     If the plant is past plant-death age, then plant-death-induced
                     mortality replaces normal growth and density dependence.
                    */
                    if (plant_ages[i] <= plant_death_age) {
                        N_out[i][t+1,j] *= exp(R[j] * (1 - Z[i]) +
                            normal_rng(0, 1) * process_error);
                    } else {
                        int after_death = plant_ages[i] - plant_death_age;
                        N_out[i][t+1,j] *= exp(log_morts[after_death] +
                            normal_rng(0, 1) * process_error);
                    }

                    // Add dispersal:
                    N_out[i][t+1,j] += (immigration - emigration);

                    // Check for extinction:
                    if (N_out[i][t+1,j] < 1) {
                        N_out[i][t+1,j] = 0;
                        extinct[i,j] = 1;
                    }

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
                        aphid_pool_j += N_out[i][t+1,j];
                        N_out[i][t+1,j] = 0;
                        extinct[i,j] = 1;
                    }
                    // Now add a portion to non-replaced plant(s):
                    aphid_pool_j /= n_non_repl_plants_;
                    for (i_ in 1:n_non_repl_plants) {
                        int i = non_repl_plants[i_];
                        // Check that it's not adding < 1 aphid to an extinct plant:
                        if (extinct[i,j] == 1 && aphid_pool_j < 1) continue;
                        // If not, add `aphid_pool_j`:
                        N_out[i][t+1,j] += aphid_pool_j;
                    }
                }
                // Now adjust index for `repl_times`
                repl_ind += 1;
            }


        }

    }

}


