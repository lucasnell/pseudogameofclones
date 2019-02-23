#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4full_model_mod) {


    class_<rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> >("model_full_model")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_full_model_namespace::model_full_model, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4full_model_plant_death_mod) {


    class_<rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> >("model_full_model_plant_death")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_full_model_plant_death_namespace::model_full_model_plant_death, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4no_among_alpha_mod) {


    class_<rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> >("model_no_among_alpha")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_no_among_alpha_namespace::model_no_among_alpha, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4no_within_alpha_mod) {


    class_<rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> >("model_no_within_alpha")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_no_within_alpha_namespace::model_no_within_alpha, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4one_alpha_mod) {


    class_<rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> >("model_one_alpha")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_one_alpha_namespace::model_one_alpha, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4one_r_mod) {


    class_<rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> >("model_one_r")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_one_r_namespace::model_one_r, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4one_r_alpha_mod) {


    class_<rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> >("model_one_r_alpha")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_one_r_alpha_namespace::model_one_r_alpha, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4pass_sigma_epsilon_mod) {


    class_<rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> >("model_pass_sigma_epsilon")

    .constructor<SEXP,SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_pass_sigma_epsilon_namespace::model_pass_sigma_epsilon, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
