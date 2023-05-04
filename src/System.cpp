
#include "System.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
void System::load_data(Rcpp::List args_data) {
  
  if (!silent) {
    print("Loading data");
  }
  
  // unpack drug concentration data
  Rcpp::List data_drug = args_data["data_drug"];
  vector<int> ind = rcpp_to_vector_int(data_drug["individual"]);
  vector<int> time = rcpp_to_vector_int(data_drug["time"]);
  vector<double> drug_conc_raw = rcpp_to_vector_double(data_drug["drug_value.conc"]);
  
  // unpack data on individual group weights
  ind_weight = rcpp_to_vector_double(args_data["ind_weight"]);
  n_ind = ind_weight.size();
  n_time = drug_conc_raw.size() / n_ind;
  
  // unpack EIR scaling factors
  eir_adjustment = rcpp_to_vector_double(args_data["eir_adjustment"]);
  eir_unique = rcpp_to_vector_double(args_data["eir_unique"]);
  eir_weight = rcpp_to_vector_double(args_data["eir_weight"]);
  
  // unpack trial data
  Rcpp::List data_control = args_data["data_control"];
  control_n = rcpp_to_vector_int(data_control["n_patients"]);
  control_n_inf = rcpp_to_vector_int(data_control["n_infected"]);
  control_time0 = rcpp_to_vector_int(data_control["time"]);
  control_time1 = rcpp_to_vector_int(data_control["time.1"]);
  
  Rcpp::List data_treat = args_data["data_treat"];
  treat_n = rcpp_to_vector_int(data_treat["n_patients"]);
  treat_n_inf = rcpp_to_vector_int(data_treat["n_infected"]);
  treat_time0 = rcpp_to_vector_int(data_treat["time"]);
  treat_time1 = rcpp_to_vector_int(data_treat["time.1"]);
  
  // reformat drug concentration data into matrix
  drug_conc = vector<vector<double>>(n_ind, vector<double>(n_time));
  for (size_t i = 0; i < ind.size(); ++i) {
    drug_conc[ind[i] - 1][time[i]] = drug_conc_raw[i];
  }
  
  // lambda weekly weighting
  n_weeks = rcpp_to_int(args_data["n_weeks"]);
  control_lambda_index = rcpp_to_matrix_int(args_data["control_lambda_index"]);
  control_lambda_weight = rcpp_to_matrix_double(args_data["control_lambda_weight"]);
}

//------------------------------------------------
void System::load_params(Rcpp::List args_params) {
  
  // MCMC parameters
  burnin = rcpp_to_int(args_params["burnin"]);
  samples = rcpp_to_int(args_params["samples"]);
  rungs = rcpp_to_int(args_params["rungs"]);
  coupling_on = rcpp_to_bool(args_params["coupling_on"]);
  beta = rcpp_to_vector_double(args_params["beta"]);
  chain = rcpp_to_int(args_params["chain"]);
  target_acceptance = rcpp_to_double(args_params["target_acceptance"]);
  
  // misc parameters
  pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  silent = rcpp_to_bool(args_params["silent"]);
}

