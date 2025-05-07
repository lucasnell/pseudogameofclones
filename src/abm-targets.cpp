
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <pcg/pcg_random.hpp>   // pcg prng

#ifndef RCPPTHREAD_OVERRIDE_COUT
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#endif
#ifndef RCPPTHREAD_OVERRIDE_CERR
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
#endif
// #ifndef RCPPTHREAD_OVERRIDE_THREAD
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override
// #endif
#include <RcppThread.h>         // multithreading


#include "pseudogameofclones_types.hpp"  // integer types
#include "pcg.hpp"              // runif_01, seed_rng functions
#include "abm-targets.hpp"      // DimensionConverter and LocationSampler classes
#include "util.hpp"      // thread_check



using namespace Rcpp;



//' Simulate targets of different types
//'
//' Simulate locations of targets of different types along an evenly spaced
//' grid of integers, where the placement of one target can affect subsequent
//' placement of other targets.
//' Target locations are drawn from all combinations of `1` to
//' `x_size-1` and `1` to `y_size-1`.
//' I'm subtracting 1 because the landscape in `searcher_sims` is defined by
//' the bounds `c(0, x_size)` and `c(0, y_size)`, so this subtraction
//' keeps a gap of 1 between the most extreme locations and the bounds
//' of the landscape.
//' This prevents targets from being located on the bounds.
//'
//'
//' @param x_size Single integer indicating x dimension of search area.
//'     Locations will be drawn from `1` to `x_size-1`.
//'     See description above for why.
//' @param y_size Single integer indicating y dimension of search area.
//'     Locations will be drawn from `1` to `y_size-1`.
//'     See description above for why.
//' @param wt_mat Square numeric matrix indicating how sample weighting on
//'     neighboring locations is affected by a target of each type being
//'     placed in a particular spot.
//'     Item `wt_mat[i,j]` indicates the effect of target type `i` on
//'     subsequent samplings of target type `j`.
//'     Values above 1 cause neighboring locations to be more likely to be
//'     sampled later, while values below 1 cause them to be less likely
//'     sampled.
//'     For locations that have been adjusted using `wt_mat` multiple times,
//'     the weights are multiplied by each other
//'     (e.g., `w *= wt_mat[1,3]` at time `t`, then
//'     `w *= wt_mat[1,2]` at time `t+1`).
//'     All weights start with values of `1`.
//'     The matrix should have the same number of rows and columns as the
//'     number of target types.
//' @param n_samples Integer vector indicating the number of samples to
//'     produce per target type. The length should equal the number of
//'     target types, and all values should be `>=1`.
//' @param allow_overlap Single logical indicating whether to allow
//'     targets to be more than one type.
//'     An example of a target being multiple types would be a plant that
//'     hosts an epiphytic bacteria and is also infected with a virus.
//'     Defaults to `TRUE`.
//' @param fill_all Single logical indicating whether to have a target at all
//'     points in the landscape. If `TRUE`, then locations where no targets
//'     were simulated will be assigned to a target type `length(n_samples)+1`.
//'     Defaults to `TRUE`.
//' @param n_lands Single integer indicating the number of independent
//'     landscapes to simulate. Each landscape will have a number of samples
//'     per type according to `n_samples`.
//'     Must be `> 0` and `< 1e6`.
//'     Defaults to `1`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Ignored if `n_lands == 1`.
//'     Defaults to `1L`.
//'
//' @return A list of length `n_lands`, where each item is a
//'     [`tibble`][tibble::tbl_df] with columns `x`, `y`, and `type`.
//'
//'
//' @export
//'
//[[Rcpp::export]]
List target_type_sims(int x_size,
                      int y_size,
                      const arma::mat& wt_mat,
                      const arma::ivec& n_samples,
                      const bool& allow_overlap = true,
                      const bool& fill_all = true,
                      const uint32& n_lands = 1,
                      const bool& show_progress = false,
                      uint32 n_threads = 1) {

    if (x_size < 2) stop("x_size must be >= 2");
    if (y_size < 2) stop("y_size must be >= 2");
    if (wt_mat.n_elem == 0) stop("wt_mat cannot be empty");
    if (!wt_mat.is_square()) stop("wt_mat must be square");
    if (arma::any(arma::vectorise(wt_mat) < 0)) stop("wt_mat cannot contain values < 0");
    if (n_samples.n_elem != wt_mat.n_rows)
        stop("length(n_samples) == nrow(wt_mat) must be true");
    if (arma::any(n_samples <= 0)) stop("n_samples cannot contain values <= 0");
    if (n_lands < 1) stop("n_lands must be > 0");
    if (n_lands >= 1e6) stop("n_lands must be < 1e6");
    thread_check(n_threads); // Check that # threads isn't too high

    // I'm reducing these by 1 bc I want to sample from 1 to floor(x_size)
    // and floor(y_size). See description in docs above for why.
    x_size--;
    y_size--;

    if (arma::any(n_samples > x_size * y_size))
        stop("n_samples cannot contain values > (x_size-1) * (y_size-1)");
    if (!allow_overlap && arma::accu(n_samples) > (x_size * y_size)) {
        stop("sum(n_samples) cannot be > (x_size-1) * (y_size-1) when allow_overlap = FALSE");
    }


    std::vector<LandSimmer> simmers;
    simmers.reserve(n_lands);
    for (uint32 i = 0; i < n_lands; i++) {
        simmers.push_back(LandSimmer(wt_mat, n_samples, x_size, y_size,
                                     allow_overlap));
    }

    RcppThread::ProgressBar prog_bar(n_lands * arma::accu(n_samples), 1);

    if (n_threads > 1U && n_lands > 1U) {
        RcppThread::parallelFor(0, n_lands, [&] (uint32 i) {
            simmers[i].run(prog_bar, show_progress);
        }, n_threads);
    } else {
        for (uint32 i = 0; i < n_lands; i++) {
            simmers[i].run(prog_bar, show_progress);
        }
    }

    List out(n_lands);
    DataFrame out_df;
    for (uint32 i = 0; i < n_lands; i++) {
        out_df = simmers[i].create_output(fill_all);
        out[i] = out_df;
    }

    return out;

}
