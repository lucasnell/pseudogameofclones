
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
