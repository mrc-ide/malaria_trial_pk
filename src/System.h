
#pragma once

#include <Rcpp.h>

#include <vector>

//------------------------------------------------
// class holding all data, parameters and functions
class System {
  
public:
  // PUBLIC OBJECTS
  
  // data
  std::vector<std::vector<double>> drug_conc;
  std::vector<double> ind_weight;
  int n_ind;
  int n_time;
  
  std::vector<int> control_n;
  std::vector<int> control_n_inf;
  std::vector<int> control_time0;
  std::vector<int> control_time1;
  
  std::vector<int> treat_n;
  std::vector<int> treat_n_inf;
  std::vector<int> treat_time0;
  std::vector<int> treat_time1;
  
  // lambda weekly weighting
  int n_weeks;
  std::vector<std::vector<int>> control_lambda_index;
  std::vector<std::vector<double>> control_lambda_weight;
  
  // MCMC parameters
  int burnin;
  int samples;
  int rungs;
  bool coupling_on;
  std::vector<double> beta;
  int chain;
  double target_acceptance;
  
  // misc parameters
  bool pb_markdown;
  bool silent;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args);
};
