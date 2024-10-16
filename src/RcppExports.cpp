// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "pseudogameofclones_types.hpp"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// logit
NumericVector logit(NumericVector p);
RcppExport SEXP _pseudogameofclones_logit(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logit(p));
    return rcpp_result_gen;
END_RCPP
}
// inv_logit
NumericVector inv_logit(NumericVector a);
RcppExport SEXP _pseudogameofclones_inv_logit(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_logit(a));
    return rcpp_result_gen;
END_RCPP
}
// leslie_matrix
NumericMatrix leslie_matrix(IntegerVector instar_days, const double& surv_juv, NumericVector surv_adult, NumericVector repro);
RcppExport SEXP _pseudogameofclones_leslie_matrix(SEXP instar_daysSEXP, SEXP surv_juvSEXP, SEXP surv_adultSEXP, SEXP reproSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type instar_days(instar_daysSEXP);
    Rcpp::traits::input_parameter< const double& >::type surv_juv(surv_juvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type surv_adult(surv_adultSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type repro(reproSEXP);
    rcpp_result_gen = Rcpp::wrap(leslie_matrix(instar_days, surv_juv, surv_adult, repro));
    return rcpp_result_gen;
END_RCPP
}
// sad_leslie
NumericVector sad_leslie(NumericMatrix leslie);
RcppExport SEXP _pseudogameofclones_sad_leslie(SEXP leslieSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type leslie(leslieSEXP);
    rcpp_result_gen = Rcpp::wrap(sad_leslie(leslie));
    return rcpp_result_gen;
END_RCPP
}
// fields_to_data_frames
List fields_to_data_frames(SEXP all_fields_ptr);
RcppExport SEXP _pseudogameofclones_fields_to_data_frames(SEXP all_fields_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type all_fields_ptr(all_fields_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(fields_to_data_frames(all_fields_ptr));
    return rcpp_result_gen;
END_RCPP
}
// fields_from_vectors
void fields_from_vectors(SEXP all_fields_ptr, std::vector<std::vector<double>>& N_vecs);
RcppExport SEXP _pseudogameofclones_fields_from_vectors(SEXP all_fields_ptrSEXP, SEXP N_vecsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type all_fields_ptr(all_fields_ptrSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<double>>& >::type N_vecs(N_vecsSEXP);
    fields_from_vectors(all_fields_ptr, N_vecs);
    return R_NilValue;
END_RCPP
}
// using_openmp
bool using_openmp();
RcppExport SEXP _pseudogameofclones_using_openmp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(using_openmp());
    return rcpp_result_gen;
END_RCPP
}
// sim_pseudogameofclones_cpp
List sim_pseudogameofclones_cpp(const uint32& n_reps, const uint32& n_fields, const std::deque<uint32>& check_for_clear, const double& clear_surv, const uint32& max_t, const uint32& save_every, const double& mean_K, const double& sd_K, const std::vector<double>& K_y_mult, const double& wilted_N, const double& wilted_mort, const arma::mat& attack_surv, const arma::mat& attack_mumm, const bool& aphid_demog_error, const bool& wasp_demog_error, const double& sigma_x, const double& sigma_y, const double& rho, const double& extinct_N, const std::vector<std::string>& aphid_name, const std::vector<arma::cube>& leslie_mat, const std::vector<arma::cube>& aphid_density_0, const std::vector<double>& alate_b0, const std::vector<double>& alate_b1, const double& alate_field_disp_p, const std::vector<double>& aphid_plant_disp_p, const std::vector<double>& plant_disp_mort, const std::vector<uint32>& field_disp_start, const std::vector<uint32>& plant_disp_start, const std::vector<uint32>& living_days, const std::vector<double>& pred_rate, const arma::mat& mum_density_0, const double& mum_smooth, const double& max_mum_density, const arma::vec& rel_attack, const double& a, const double& k, const double& h, const double& wasp_badger_n, const std::vector<double>& wasp_density_0, const std::vector<uint32>& wasp_delay, const double& wasp_disp_m0, const double& wasp_disp_m1, const std::vector<double>& wasp_field_attract, const double& sex_ratio, const std::vector<double>& s_y, const std::vector<bool>& constant_wasps, const std::vector<uint32>& perturb_when, const std::vector<uint32>& perturb_where, const std::vector<uint32>& perturb_who, const std::vector<double>& perturb_how, const arma::umat& extra_plant_removals_mat, const bool& sep_adults, uint32 n_threads, const bool& show_progress);
RcppExport SEXP _pseudogameofclones_sim_pseudogameofclones_cpp(SEXP n_repsSEXP, SEXP n_fieldsSEXP, SEXP check_for_clearSEXP, SEXP clear_survSEXP, SEXP max_tSEXP, SEXP save_everySEXP, SEXP mean_KSEXP, SEXP sd_KSEXP, SEXP K_y_multSEXP, SEXP wilted_NSEXP, SEXP wilted_mortSEXP, SEXP attack_survSEXP, SEXP attack_mummSEXP, SEXP aphid_demog_errorSEXP, SEXP wasp_demog_errorSEXP, SEXP sigma_xSEXP, SEXP sigma_ySEXP, SEXP rhoSEXP, SEXP extinct_NSEXP, SEXP aphid_nameSEXP, SEXP leslie_matSEXP, SEXP aphid_density_0SEXP, SEXP alate_b0SEXP, SEXP alate_b1SEXP, SEXP alate_field_disp_pSEXP, SEXP aphid_plant_disp_pSEXP, SEXP plant_disp_mortSEXP, SEXP field_disp_startSEXP, SEXP plant_disp_startSEXP, SEXP living_daysSEXP, SEXP pred_rateSEXP, SEXP mum_density_0SEXP, SEXP mum_smoothSEXP, SEXP max_mum_densitySEXP, SEXP rel_attackSEXP, SEXP aSEXP, SEXP kSEXP, SEXP hSEXP, SEXP wasp_badger_nSEXP, SEXP wasp_density_0SEXP, SEXP wasp_delaySEXP, SEXP wasp_disp_m0SEXP, SEXP wasp_disp_m1SEXP, SEXP wasp_field_attractSEXP, SEXP sex_ratioSEXP, SEXP s_ySEXP, SEXP constant_waspsSEXP, SEXP perturb_whenSEXP, SEXP perturb_whereSEXP, SEXP perturb_whoSEXP, SEXP perturb_howSEXP, SEXP extra_plant_removals_matSEXP, SEXP sep_adultsSEXP, SEXP n_threadsSEXP, SEXP show_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const uint32& >::type n_reps(n_repsSEXP);
    Rcpp::traits::input_parameter< const uint32& >::type n_fields(n_fieldsSEXP);
    Rcpp::traits::input_parameter< const std::deque<uint32>& >::type check_for_clear(check_for_clearSEXP);
    Rcpp::traits::input_parameter< const double& >::type clear_surv(clear_survSEXP);
    Rcpp::traits::input_parameter< const uint32& >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< const uint32& >::type save_every(save_everySEXP);
    Rcpp::traits::input_parameter< const double& >::type mean_K(mean_KSEXP);
    Rcpp::traits::input_parameter< const double& >::type sd_K(sd_KSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type K_y_mult(K_y_multSEXP);
    Rcpp::traits::input_parameter< const double& >::type wilted_N(wilted_NSEXP);
    Rcpp::traits::input_parameter< const double& >::type wilted_mort(wilted_mortSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type attack_surv(attack_survSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type attack_mumm(attack_mummSEXP);
    Rcpp::traits::input_parameter< const bool& >::type aphid_demog_error(aphid_demog_errorSEXP);
    Rcpp::traits::input_parameter< const bool& >::type wasp_demog_error(wasp_demog_errorSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_x(sigma_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma_y(sigma_ySEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double& >::type extinct_N(extinct_NSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type aphid_name(aphid_nameSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::cube>& >::type leslie_mat(leslie_matSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::cube>& >::type aphid_density_0(aphid_density_0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alate_b0(alate_b0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alate_b1(alate_b1SEXP);
    Rcpp::traits::input_parameter< const double& >::type alate_field_disp_p(alate_field_disp_pSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type aphid_plant_disp_p(aphid_plant_disp_pSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type plant_disp_mort(plant_disp_mortSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type field_disp_start(field_disp_startSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type plant_disp_start(plant_disp_startSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type living_days(living_daysSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type pred_rate(pred_rateSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mum_density_0(mum_density_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type mum_smooth(mum_smoothSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_mum_density(max_mum_densitySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rel_attack(rel_attackSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const double& >::type wasp_badger_n(wasp_badger_nSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type wasp_density_0(wasp_density_0SEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type wasp_delay(wasp_delaySEXP);
    Rcpp::traits::input_parameter< const double& >::type wasp_disp_m0(wasp_disp_m0SEXP);
    Rcpp::traits::input_parameter< const double& >::type wasp_disp_m1(wasp_disp_m1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type wasp_field_attract(wasp_field_attractSEXP);
    Rcpp::traits::input_parameter< const double& >::type sex_ratio(sex_ratioSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type s_y(s_ySEXP);
    Rcpp::traits::input_parameter< const std::vector<bool>& >::type constant_wasps(constant_waspsSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type perturb_when(perturb_whenSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type perturb_where(perturb_whereSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type perturb_who(perturb_whoSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type perturb_how(perturb_howSEXP);
    Rcpp::traits::input_parameter< const arma::umat& >::type extra_plant_removals_mat(extra_plant_removals_matSEXP);
    Rcpp::traits::input_parameter< const bool& >::type sep_adults(sep_adultsSEXP);
    Rcpp::traits::input_parameter< uint32 >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show_progress(show_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_pseudogameofclones_cpp(n_reps, n_fields, check_for_clear, clear_surv, max_t, save_every, mean_K, sd_K, K_y_mult, wilted_N, wilted_mort, attack_surv, attack_mumm, aphid_demog_error, wasp_demog_error, sigma_x, sigma_y, rho, extinct_N, aphid_name, leslie_mat, aphid_density_0, alate_b0, alate_b1, alate_field_disp_p, aphid_plant_disp_p, plant_disp_mort, field_disp_start, plant_disp_start, living_days, pred_rate, mum_density_0, mum_smooth, max_mum_density, rel_attack, a, k, h, wasp_badger_n, wasp_density_0, wasp_delay, wasp_disp_m0, wasp_disp_m1, wasp_field_attract, sex_ratio, s_y, constant_wasps, perturb_when, perturb_where, perturb_who, perturb_how, extra_plant_removals_mat, sep_adults, n_threads, show_progress));
    return rcpp_result_gen;
END_RCPP
}
// restart_fill_other_pars
SEXP restart_fill_other_pars(SEXP all_fields_in_ptr, const double& K, const std::vector<double>& alate_b0, const std::vector<double>& alate_b1, const double& alate_field_disp_p, const std::vector<double>& K_y_mult, const std::vector<double>& s_y, const double& a, const double& k, const double& h, const double& wasp_disp_m0, const double& wasp_disp_m1, const std::vector<double>& wasp_field_attract, const double& mum_smooth, const std::vector<double>& pred_rate);
RcppExport SEXP _pseudogameofclones_restart_fill_other_pars(SEXP all_fields_in_ptrSEXP, SEXP KSEXP, SEXP alate_b0SEXP, SEXP alate_b1SEXP, SEXP alate_field_disp_pSEXP, SEXP K_y_multSEXP, SEXP s_ySEXP, SEXP aSEXP, SEXP kSEXP, SEXP hSEXP, SEXP wasp_disp_m0SEXP, SEXP wasp_disp_m1SEXP, SEXP wasp_field_attractSEXP, SEXP mum_smoothSEXP, SEXP pred_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type all_fields_in_ptr(all_fields_in_ptrSEXP);
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alate_b0(alate_b0SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alate_b1(alate_b1SEXP);
    Rcpp::traits::input_parameter< const double& >::type alate_field_disp_p(alate_field_disp_pSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type K_y_mult(K_y_multSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type s_y(s_ySEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type h(hSEXP);
    Rcpp::traits::input_parameter< const double& >::type wasp_disp_m0(wasp_disp_m0SEXP);
    Rcpp::traits::input_parameter< const double& >::type wasp_disp_m1(wasp_disp_m1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type wasp_field_attract(wasp_field_attractSEXP);
    Rcpp::traits::input_parameter< const double& >::type mum_smooth(mum_smoothSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type pred_rate(pred_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(restart_fill_other_pars(all_fields_in_ptr, K, alate_b0, alate_b1, alate_field_disp_p, K_y_mult, s_y, a, k, h, wasp_disp_m0, wasp_disp_m1, wasp_field_attract, mum_smooth, pred_rate));
    return rcpp_result_gen;
END_RCPP
}
// restart_experiments_cpp
List restart_experiments_cpp(SEXP all_fields_ptr, const uint32& max_t, const uint32& save_every, const std::deque<uint32>& check_for_clear, const bool& stage_ts_out, const bool& sep_adults, const bool& show_progress, const std::vector<uint32>& perturb_when, const std::vector<uint32>& perturb_where, const std::vector<uint32>& perturb_who, const std::vector<double>& perturb_how);
RcppExport SEXP _pseudogameofclones_restart_experiments_cpp(SEXP all_fields_ptrSEXP, SEXP max_tSEXP, SEXP save_everySEXP, SEXP check_for_clearSEXP, SEXP stage_ts_outSEXP, SEXP sep_adultsSEXP, SEXP show_progressSEXP, SEXP perturb_whenSEXP, SEXP perturb_whereSEXP, SEXP perturb_whoSEXP, SEXP perturb_howSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type all_fields_ptr(all_fields_ptrSEXP);
    Rcpp::traits::input_parameter< const uint32& >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< const uint32& >::type save_every(save_everySEXP);
    Rcpp::traits::input_parameter< const std::deque<uint32>& >::type check_for_clear(check_for_clearSEXP);
    Rcpp::traits::input_parameter< const bool& >::type stage_ts_out(stage_ts_outSEXP);
    Rcpp::traits::input_parameter< const bool& >::type sep_adults(sep_adultsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type show_progress(show_progressSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type perturb_when(perturb_whenSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type perturb_where(perturb_whereSEXP);
    Rcpp::traits::input_parameter< const std::vector<uint32>& >::type perturb_who(perturb_whoSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type perturb_how(perturb_howSEXP);
    rcpp_result_gen = Rcpp::wrap(restart_experiments_cpp(all_fields_ptr, max_t, save_every, check_for_clear, stage_ts_out, sep_adults, show_progress, perturb_when, perturb_where, perturb_who, perturb_how));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pseudogameofclones_logit", (DL_FUNC) &_pseudogameofclones_logit, 1},
    {"_pseudogameofclones_inv_logit", (DL_FUNC) &_pseudogameofclones_inv_logit, 1},
    {"_pseudogameofclones_leslie_matrix", (DL_FUNC) &_pseudogameofclones_leslie_matrix, 4},
    {"_pseudogameofclones_sad_leslie", (DL_FUNC) &_pseudogameofclones_sad_leslie, 1},
    {"_pseudogameofclones_fields_to_data_frames", (DL_FUNC) &_pseudogameofclones_fields_to_data_frames, 1},
    {"_pseudogameofclones_fields_from_vectors", (DL_FUNC) &_pseudogameofclones_fields_from_vectors, 2},
    {"_pseudogameofclones_using_openmp", (DL_FUNC) &_pseudogameofclones_using_openmp, 0},
    {"_pseudogameofclones_sim_pseudogameofclones_cpp", (DL_FUNC) &_pseudogameofclones_sim_pseudogameofclones_cpp, 55},
    {"_pseudogameofclones_restart_fill_other_pars", (DL_FUNC) &_pseudogameofclones_restart_fill_other_pars, 15},
    {"_pseudogameofclones_restart_experiments_cpp", (DL_FUNC) &_pseudogameofclones_restart_experiments_cpp, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_pseudogameofclones(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
