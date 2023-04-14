
#pragma once

#include "System.h"
#include "misc.h"

#include <Rcpp.h>

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  System * s_ptr;
  
  // thermodynamic power
  double beta;
  
  // parameters
  std::vector<double> lambda;
  double min_prob;
  double half_point;
  double hill_power;
  
  // intermediate objects
  std::vector<std::vector<double>> drug_pow;
  std::vector<std::vector<double>> drug_pow_prop;
  std::vector<std::vector<double>> exp_rate;
  
  // proposal parameters
  double bw_lambda;
  double bw_min_prob;
  double bw_half_point;
  double bw_hill_power;
  double bw_stepsize;
  
  // likelihoods and priors
  double loglike_control;
  double loglike_treat;
  double loglike;
  double logprior;
  
  // store acceptance rates
  int accept_count;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialise everything except for likelihood and prior values
  void init(System &s, double beta);
  
  // main methods
  void update(bool RM_on, int iteration);
  void update_lambda(bool RM_on, int iteration);
  void update_min_prob(bool RM_on, int iteration);
  void update_half_point(bool RM_on, int iteration);
  void update_hill_power(bool RM_on, int iteration);
  void recalc_drug_pow(std::vector<std::vector<double>> &mat, double k);
  double get_loglike_control(std::vector<double> &lambda_);
  double get_loglike_treat(std::vector<double> &lambda_,
                           double min_prob_, double half_point_,
                           std::vector<std::vector<double>> &drug_pow_);
  double get_logprior(std::vector<double> &lambda_, double min_prob_,
                      double half_point_, double hill_power_);
  
};
