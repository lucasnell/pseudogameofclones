// I no longer want this to compile
# ifdef __CLONEWARS_PLANTS_H
// # ifndef __CLONEWARS_PLANTS_H
# define __CLONEWARS_PLANTS_H


#include <RcppArmadillo.h>
#include <random> // C++11 distributions
#include <pcg/pcg_random.hpp> // pcg PRNG
#include <progress.hpp>  // for the progress bar

#ifdef _OPENMP
#include <omp.h>
#endif

#include "clonewars_types.hpp"
#include "pcg.hpp"
#include "simulate.hpp"

using namespace Rcpp;





/*
 ========================================================================================
 ========================================================================================

 Functions for OnePatch class

 ========================================================================================
 ========================================================================================
 */


/*
 Reset back to original values.
 */
void OnePatch::reset(const double& zeta_) {
    Nt = N0;
    Nt1 = N0;
    extinct = arma::Row<unsigned short>(N0.n_elem, arma::fill::zeros);
    zeta = zeta_;
    DD = arma::as_scalar(A * Nt.t()) / mean_alpha;
    Z = arma::sum(Nt);
    age = 0;
    empty = arma::accu(N0) == 0;
    return;
}

/*
 Clear to no aphids and a new zeta
*/
void OnePatch::clear(const double& zeta_) {
    Nt.fill(0);
    Nt1.fill(0);
    extinct.fill(1);
    zeta = zeta_;
    DD = 0;
    Z = 0;
    age = 0;
    empty = true;
    return;
}

/*
 Emigration of one line from this patch to all other patches:
 */
arma::vec OnePatch::emigration(const uint32& j,
                               const arma::mat& D_vec,
                               pcg32& eng) {

    arma::vec disp_out(n_patches, arma::fill::zeros);

    const double& Ntj(Nt(j));

    if (Ntj == 0 || n_patches == 1) return disp_out;

    // Else we have to sample from Poisson:
    // Calculate lambda:
    double lambda_ = D_vec(j) * Ntj;
    // Because we're splitting mu_ over the number of patches:
    double np = static_cast<double>(n_patches);
    lambda_ /= np;
    // Update Poisson distribution based on this:
    poisson_distr.param(
        std::poisson_distribution<arma::sword>::param_type(lambda_));
    for (uint32 ii = 0; ii < n_patches; ii++) {
        if (ii == i) {
            /*
            No point in sampling this combination bc it's dispersal back to
            the same patch, so summing immigrants and emigrants would cancel
            out anyway.
            */
            ;
        } else {
            disp_out(ii) = static_cast<double>(poisson_distr(eng));
        }
    }
    // Making absolutely sure that dispersal never exceeds the number possible:
    double total_emigrants = arma::accu(disp_out);
    if (total_emigrants > Ntj) {
        double extras = total_emigrants - Ntj;
        std::vector<uint32> extra_inds;
        extra_inds.reserve(n_patches);
        for (uint32 ii = 0; ii < n_patches; ii++) {
            if (disp_out(ii) > 0) extra_inds.push_back(ii);
        }
        while (extras > 0) {
            uint32 rnd = runif_01(eng) * extra_inds.size();
            disp_out(extra_inds[rnd])--;
            if (disp_out(extra_inds[rnd]) == 0) {
                extra_inds.erase(extra_inds.begin() + rnd);
            }
            extras--;
        }
    }

    return disp_out;
}

/*
Same thing as above, but overloaded for not including dispersal stochasticity
*/
arma::vec OnePatch::emigration(const uint32& j, const arma::mat& D_vec) {

    const double& Ntj(Nt(j));

    if (Ntj == 0 || n_patches == 1) return arma::vec(n_patches, arma::fill::zeros);

    // Else we have to calculate E(# dispersed)
    // Calculate lambda:
    double lambda_ = D_vec(j) * Ntj;
    // Because we're splitting mu_ over the number of patches:
    double np = static_cast<double>(Nt.n_rows);
    lambda_ /= np;

    arma::vec disp_out(n_patches);
    disp_out.fill(lambda_);
    disp_out(i) = 0;

    return disp_out;
}




// Check to see if this patch will be replaced for a given threshold
void OnePatch::replace_check(const double& repl_threshold,
                             const double& zeta_t_thresh,
                             std::vector<uint32>& repl_inds,
                             double& not_replaced_Z) {
    // if (Z > repl_threshold || (empty && age > 0) || (zeta_t > zeta_t_thresh)) {
    if (Z > repl_threshold || (zeta_t > zeta_t_thresh)) {
        repl_inds.push_back(this->i);
    } else not_replaced_Z += Z;
    return;
}


/*
 Iterate one time step
*/
void OnePatch::iterate(const arma::rowvec& R,
                       const arma::mat& emigrants,
                       const arma::mat& immigrants,
                       const double& extinct_N,
                       const double& mu_time,
                       pcg32& eng) {

    // Update density dependency:
    DD = arma::as_scalar(A * Nt.t()) / mean_alpha;
    // Overall effect of patch health deterioration over time:
    zeta_t = std::exp(zeta * (age - mu_time));

    for (uint32 j = 0; j < n_lines; j++) {

        // Dispersal of this aphid line to and from other patches.
        double imm = immigrants(i, j);
        double em = emigrants(i, j);

        // If it's extinct and no one's coming in, make 0 and skip the rest:
        if (extinct(j) == 1 && (imm - em) <= 0.0) {
            Nt1(j) = 0;
            continue;
        }
        // If it's extinct and immigrants are coming in, adjust extinct matrix:
        if (extinct(j) == 1 && (imm - em) > 0.0) {
            extinct(j) = 0;
        }

        /*
        Add growth, density dependence, and process error.
        */
        Nt1(j) *= std::exp(R(j) * (1 - A(j) * DD * zeta_t) +
            normal_distr(eng) * process_error);

        // Add dispersal:
        Nt1(j) += (imm - em);

        // Check for extinction:
        if (Nt1(j) < extinct_N) {
            Nt1(j) = 0;
            extinct(j) = 1;
        }

    }

    Nt = Nt1;
    Z = arma::sum(Nt1);
    if (Z > 0) {
        age++;
        empty = false;
    } else empty = true;

    return;
}




// Fill part of a matrix with info from this patch:
void OnePatch::fill_matrix(arma::mat& matrix, const uint32& row_start) {

    if (matrix.n_cols == 5) {
        for (uint32 ii = 0; ii < n_lines; ii++) {
            matrix(row_start + ii, 2) = i;
            matrix(row_start + ii, 3) = ii;
            matrix(row_start + ii, 4) = Nt1(ii);
        }
    } else if (matrix.n_cols == 4) {
        matrix(row_start, 2) = i;
        matrix(row_start, 3) = Z;
    }

    return;
}






/*
 ========================================================================================
 ========================================================================================

        Functions for AllPatches class

 ========================================================================================
 ========================================================================================
 */



// Reset for another rep
void AllPatches::reset(const uint32& rep_, pcg32& eng) {
    rep = static_cast<double>(rep_);
    repl_times = repl_times0;
    uint32 chunk_size = n_patches;
    if (!by_patch) chunk_size *= n_lines;
    chunk_size *= n_saves;
    output_row = rep_ * chunk_size;
    for (OnePatch& patch : patches) patch.reset(zeta_distr.sample(eng));
    return;
}



/*
 Update emigrants and immigrants matrices
 */
void AllPatches::update_dispersal(pcg32& eng) {

    immigrants.fill(0);
    uint32 n_lines = emigrants.n_cols;
    uint32 n_patches = emigrants.n_rows;

    arma::vec e_ij;
    for (uint32 j = 0; j < n_lines; j++) {
        for (uint32 i = 0; i < n_patches; i++) {
            // Emigrants of line j from patch i to all others:
            if (disp_error) {
                e_ij = patches[i].emigration(j, D_vec, eng);
            } else {
                e_ij = patches[i].emigration(j, D_vec);
            }
            emigrants(i, j) = arma::accu(e_ij);
            immigrants.col(j) += e_ij;
        }
    }

    return;

}

/*
Replace patches when they exceed a threshold for number of aphids
*/
void AllPatches::replace_patches(const uint32& t, pcg32& eng) {

    // Do nothing if this isn't a day to replace patches
    if (repl_times.empty()) return;
    if ((t+1) != repl_times.front()) return;

    repl_times.pop_front();

    /*
    Find the patches to replace and those not to, and update patch_days for
    those that are replaced.
    */
    std::vector<uint32> replaced;
    replaced.reserve(n_patches);
    double not_replaced_Z = 0;  // total aphids in non-replaced patches
    for (uint32 i = 0; i < n_patches; i++) {
        patches[i].replace_check(repl_threshold, zeta_t_thresh, replaced, not_replaced_Z);
    }

    // If none to replace, then continue
    if (replaced.size() == 0) return;

    /*
     In the rare event that replacement will cause extinction, use half the indices.
     If it still causes extinction, remove indices until it does no longer:
    */
    if (not_replaced_Z == 0) {
        uint32 new_size = replaced.size() / 2;
        while (replaced.size() > new_size) {
            not_replaced_Z += patches[replaced.back()].Z;
            replaced.pop_back();
        }
    }
    // while (not_replaced_Z == 0 && !replaced.empty()) {
    //     not_replaced_Z += patches[replaced.back()].Z;
    //     replaced.pop_back();
    // }

    for (uint32& ii : replaced) {
        patches[ii].clear(zeta_distr.sample(eng));
    }

    return;
}




// Fill out matrix of output if necessary:
void AllPatches::fill_matrix(arma::mat& out_matrix,
                             const uint32& t1) {

    if (t1 % save_every != 0 && t1 != max_t) return;

    double t_(t1);

    if (!by_patch) {

        if (out_matrix.n_cols != 5) stop("ncol(matrix) should be 5 when !by_patch.");
        uint32 needed_rows = n_patches * n_lines;
        if (out_matrix.n_rows < (output_row + needed_rows)) {
            Rcout << out_matrix.n_rows << " (should be >=) " <<
                (output_row + needed_rows) << std::endl;
            stop("not enough rows.");
        }

        for (uint32 ii = 0; ii < n_patches; ii++) {
            out_matrix(arma::span(output_row, output_row+n_lines-1), 0).fill(rep);
            out_matrix(arma::span(output_row, output_row+n_lines-1), 1).fill(t_);
            patches[ii].fill_matrix(out_matrix, output_row);
            output_row += n_lines;
        }

    } else {

        if (out_matrix.n_cols != 4) stop("ncol(matrix) should be 4 when by_patch.");
        uint32 needed_rows = n_patches;
        if (out_matrix.n_rows < (output_row + needed_rows)) stop("not enough rows.");

        for (uint32 ii = 0; ii < n_patches; ii++) {
            out_matrix(output_row, 0) = rep;
            out_matrix(output_row, 1) = t_;
            patches[ii].fill_matrix(out_matrix, output_row);
            output_row++;
        }

    }

    return;
}




// Iterate all time steps, saving output on the way
void AllPatches::run_series(arma::mat& out_matrix, pcg32& eng) {

    // Get initial conditions
    fill_matrix(out_matrix, 0);

    for (uint32 t = 0; t < max_t; t++) {

        // Generate numbers of dispersed aphids:
        update_dispersal(eng);

        // This has the following error: "Mat::operator(): index out of bounds"
        // Update population abundances
        for (OnePatch& patch : patches) {
            patch.iterate(R, emigrants, immigrants, extinct_N, mu_time, eng);
        }

        // Summarize output if necessary:
        fill_matrix(out_matrix, t+1);

        // Replace patches if necessary:
        replace_patches(t, eng);

    }

    return;
}




//' Calculate the number of rows per rep.
//'
//' @noRd
//'
inline uint32 calc_rep_rows(const uint32& max_t,
                            const uint32& save_every,
                            const uint32& n_lines,
                            const uint32& n_patches,
                            const bool& by_patch) {

    uint32 n_saves = (max_t / save_every) + 1;
    if (max_t % save_every > 0) n_saves++;

    uint32 nrows_ = n_patches;
    if (!by_patch) nrows_ *= n_lines;
    nrows_ *= n_saves;

    return nrows_;
}





//' Simulate multiple reps.
//'
//'
//' @inheritParams sim_reps
//'
//' @noRd
//'
//'
//[[Rcpp::export]]
arma::mat sim_reps_(const uint32& n_reps,
                    const uint32& max_t,
                    const arma::mat& N0,
                    const arma::rowvec& R,
                    const arma::rowvec& A,
                    const arma::vec& D_vec,
                    const double& process_error,
                    const bool& disp_error,
                    const double& log_zeta_mean,
                    const double& log_zeta_sd,
                    const double& zeta_t_thresh,
                    const double& mu_time,
                    const std::deque<uint32>& repl_times,
                    const double& repl_threshold,
                    const double& extinct_N,
                    const uint32& save_every,
                    const bool& by_patch,
                    const uint32& n_cores,
                    const bool& show_progress) {


    const std::vector<std::vector<uint64> > seeds = mc_seeds(n_cores);

    const uint32 n_patches = N0.n_rows;
    const uint32 n_lines = N0.n_cols;

    // Number of rows per rep:
    uint32 rep_rows = calc_rep_rows(max_t, save_every, n_lines, n_patches, by_patch);
    uint32 rep_cols = 4;
    if (!by_patch) rep_cols++;

    std::string err_msg = "";

    if (R.n_elem != n_lines) {
        err_msg += "ERROR: R.n_elem != n_lines\n";
    }
    if (A.n_elem != n_lines) {
        err_msg += "ERROR: A.n_elem != n_lines\n";
    }
    if (D_vec.n_elem != n_lines) {
        err_msg += "ERROR: D_vec.n_elem != n_lines\n";
    }
    if (repl_times.size() == 0) {
        err_msg += "ERROR: repl_times.size() == 0\n";
    }
    if (n_patches <= 1) {
        err_msg += "ERROR: n_patches <= 1\n";
    }

    if (err_msg.size() > 0) throw(Rcpp::exception(err_msg.c_str(), false));


    if (show_progress) Rcout << "Starting simulations..." << std::endl;

    Progress p(n_reps, show_progress);


    /*
     --------------
     Create output object
     --------------
     */
    arma::mat out_matrix(n_reps * rep_rows, rep_cols);

    #ifdef _OPENMP
    #pragma omp parallel default(shared) num_threads(n_cores) if(n_cores > 1)
    {
    #endif

    std::vector<uint64> active_seeds;

    // Write the active seed per core or just write one of the seeds.
    #ifdef _OPENMP
    uint32 active_thread = omp_get_thread_num();
    active_seeds = seeds[active_thread];
    #else
    active_seeds = seeds[0];
    #endif

    // pcg prng
    pcg32 eng = seeded_pcg(active_seeds);

    AllPatches metapop(max_t, N0, R, A, D_vec, process_error, disp_error,
                       log_zeta_mean, log_zeta_sd, zeta_t_thresh, mu_time, repl_times,
                       repl_threshold, extinct_N, save_every, by_patch, 0.0, eng);

    // Parallelize the Loop
    #ifdef _OPENMP
    #pragma omp for schedule(static)
    #endif
    for (uint32 r = 0; r < n_reps; r++) {
        if (p.is_aborted()) continue;

        metapop.reset(r, eng);
        metapop.run_series(out_matrix, eng);

        p.increment();
    }

    #ifdef _OPENMP
    }
    #endif

    if (show_progress) Rcout << "... finished!" << std::endl;


    return out_matrix;

}


#endif
