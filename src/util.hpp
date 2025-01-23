#ifndef __PSEUDOGAMEOFCLONES_UTIL_H
#define __PSEUDOGAMEOFCLONES_UTIL_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <RcppThread.h>         // multithreading

#include "pseudogameofclones_types.hpp"


using namespace Rcpp;



//' Check that the number of threads doesn't exceed the number available.
 //'
 //' @noRd
 //'
 inline void thread_check(uint32& n_threads) {

     if (n_threads == 0) n_threads = 1;

     uint32 max_threads = std::thread::hardware_concurrency();

     if (n_threads > max_threads) {
         std::string mt_str = std::to_string(max_threads);
         std::string err_msg = "\nThe number of requested threads (" +
             std::to_string(n_threads) +
             ") exceeds the max available on the system (" + mt_str + ").";
         stop(err_msg.c_str());
     }

     return;
 }



#endif
