
#include "System.h"
#include "Particle.h"

#include <Rcpp.h>


//------------------------------------------------
// main Rcpp function, deployed from R
// [[Rcpp::export]]
Rcpp::List main_cpp(Rcpp::List args);

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(std::vector<Particle> &particle_vec, std::vector<int> &mc_accept);
