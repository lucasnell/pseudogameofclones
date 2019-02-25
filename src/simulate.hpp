# ifndef __CLONEWARS_SIMULATE_H
# define __CLONEWARS_SIMULATE_H


#include <RcppArmadillo.h>
#include <random> // C++11 distributions
#include <deque>  // deque class
#include <pcg/pcg_random.hpp> // pcg PRNG
#include <progress.hpp>  // for the progress bar



#include "clonewars_types.hpp"
#include "pcg.hpp"

using namespace Rcpp;



/*
 This samples for zeta values and allows for no randomness if log_zeta_sd <= 0
 */
struct ZetaSampler {

    ZetaSampler(const double& log_zeta_mean,
                const double& log_zeta_sd)
    : lnorm_distr(log_zeta_mean, log_zeta_sd),
      use_distr(log_zeta_sd > 0),
      zeta_mean(0) {

        if (!use_distr) {
            zeta_mean = std::exp(log_zeta_mean + ((log_zeta_sd * log_zeta_sd) / 2));
        }

    }

    double sample(pcg32& eng) {
        if (use_distr) return lnorm_distr(eng);
        return zeta_mean;
    }

private:

    std::lognormal_distribution<double> lnorm_distr;
    bool use_distr;
    double zeta_mean;
};






/*
One patch of continuous habitat:
*/
struct OnePatch {

    arma::rowvec Nt;                    // abundances at time t
    arma::rowvec Nt1;                   // abundances at time t+1
    arma::Row<unsigned short> extinct;  // keeping track of extinctions
    double zeta;                        // speed of patch deterioration through time
    double zeta_t;                      // current patch health: exp(zeta * (t - mu_t))
    double DD;                          // overall density dependence for this patch
    double Z;                           // overall aphid population size on this patch
    double age;                         // current patch age
    bool empty;                         // boolean for whether no aphids are on this patch

    OnePatch() {}
    OnePatch(const uint32& max_t,
             const arma::rowvec& N0_,
             const arma::rowvec& A_,
             const double& zeta_,
             const double& process_error_,
             const uint32& i_,
             const uint32& n_patches_)
        : Nt(N0_),
          Nt1(N0_),
          extinct(N0_.n_elem, arma::fill::zeros),
          zeta(zeta_),
          zeta_t(0),
          DD(0),
          Z(0),
          age(0),
          empty(arma::accu(N0_) == 0),
          poisson_distr(1.0),
          normal_distr(0, 1.0),
          process_error(process_error_),
          i(i_),
          n_patches(n_patches_),
          n_lines(N0_.n_elem),
          N0(N0_),
          A(A_),
          mean_alpha(arma::mean(A_)) {

        DD = arma::as_scalar(A * Nt.t()) / mean_alpha;
        Z = arma::sum(Nt);

    }

    /*
     Reset back to original values.
     */
    void reset(const double& zeta_);

    /*
     Clear to no aphids and a new zeta
     */
    void clear(const double& zeta_);

    /*
     Emigration of one line from this patch to all other patches:
    */
    arma::vec emigration(const uint32& j,
                         const arma::mat& D_vec,
                         pcg32& eng);

    /*
     Same thing as above, but overloaded for not including dispersal stochasticity
    */
    arma::vec emigration(const uint32& j, const arma::mat& D_vec);


    // Check to see if this patch will be replaced for a given threshold
    void replace_check(const double& repl_threshold,
                       const double& zeta_t_thresh,
                       std::vector<uint32>& repl_inds,
                       double& not_replaced_Z);

    /*
     Iterate one time step
     */
    void iterate(const arma::rowvec& R,
                 const arma::mat& emigrants,
                 const arma::mat& immigrants,
                 const double& extinct_N,
                 const double& mu_time,
                 pcg32& eng);

    // Fill part of a matrix with info from this patch:
    void fill_matrix(arma::mat& matrix, const uint32& row_start);


private:
    std::poisson_distribution<arma::sword> poisson_distr;   // for sampling emigration
    std::normal_distribution<double> normal_distr;          // for process error
    double process_error;                                   // SD of process error
    uint32 i;                                               // index for this patch
    uint32 n_patches;                                       // number of patches
    uint32 n_lines;                                         // number of lines
    arma::rowvec N0;                                        // starting N for resetting
    arma::rowvec A;                                         // density dependences (alpha)
    double mean_alpha;                                      // mean alpha among all lines

};



class AllPatches {
public:

    std::vector<OnePatch> patches;
    arma::mat emigrants;
    arma::mat immigrants;


    AllPatches(
        const uint32& max_t_,
        const arma::mat& N0_,
        const arma::rowvec& R_,
        const arma::rowvec& A_,
        const arma::vec& D_vec_,
        const double& process_error_,
        const bool& disp_error_,
        const double& log_zeta_mean_,
        const double& log_zeta_sd_,
        const double& zeta_t_thresh_,
        const double& mu_time_,
        const std::deque<uint32>& repl_times_,
        const double& repl_threshold_,
        const double& extinct_N_,
        const uint32& save_every_,
        const bool& by_patch_,
        const double& rep_,
        pcg32& eng)
    : patches(N0_.n_rows),
      emigrants(N0_.n_rows, N0_.n_cols),
      immigrants(N0_.n_rows, N0_.n_cols),
      max_t(max_t_),
      n_patches(N0_.n_rows),
      n_lines(N0_.n_cols),
      R(R_),
      A(A_),
      D_vec(D_vec_),
      process_error(process_error_),
      disp_error(disp_error_),
      mu_time(mu_time_),
      repl_times(repl_times_),
      repl_times0(repl_times_),
      repl_threshold(repl_threshold_),
      extinct_N(extinct_N_),
      save_every(save_every_),
      by_patch(by_patch_),
      rep(rep_),
      n_saves(0),
      output_row(0),
      zeta_distr(log_zeta_mean_, log_zeta_sd_),
      zeta_t_thresh(zeta_t_thresh_){

        for (uint32 i = 0; i < N0_.n_rows; i++) {
            double zeta = zeta_distr.sample(eng);
            patches[i] = OnePatch(max_t_, N0_.row(i), A, zeta, process_error,
                                  i, N0_.n_rows);
        }

        n_saves = (max_t / save_every) + 1;
        if (max_t % save_every > 0) n_saves++;

        uint32 chunk_size = n_patches;
        if (!by_patch) chunk_size *= n_lines;
        chunk_size *= n_saves;
        output_row = rep_ * chunk_size;

    }


    // Reset for another rep
    void reset(const uint32& rep_, pcg32& eng);

    /*
     Update emigrants and immigrants matrices
     */
    void update_dispersal(pcg32& eng);

    /*
     Replace patches when they exceed a threshold for number of aphids
     */
    void replace_patches(const uint32& t, pcg32& eng);

    // Fill out matrix of output if necessary:
    void fill_matrix(arma::mat& out_matrix, const uint32& t1);


    // Iterate all time steps, saving output on the way
    void run_series(arma::mat& out_matrix, pcg32& eng);



private:

    uint32 max_t;
    uint32 n_patches;
    uint32 n_lines;
    arma::rowvec R;
    arma::rowvec A;
    arma::vec D_vec;
    double process_error;
    bool disp_error;
    double mu_time;
    std::deque<uint32> repl_times;
    std::deque<uint32> repl_times0;  // original one, so that resetting this obj is easy
    double repl_threshold;
    double extinct_N;
    uint32 save_every;
    bool by_patch;
    double rep;         // rep number to use for filling output
    uint32 n_saves;     // number of time steps that will be saved
    uint32 output_row;  // keeps track of what row to start with in the output

    ZetaSampler zeta_distr;
    double zeta_t_thresh;
};






#endif
