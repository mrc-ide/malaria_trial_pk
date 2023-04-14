
#include "System.h"
#include "misc.h"

using namespace std;

void System::load(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
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
  
  if (!silent) {
    print("Loading data");
  }
  
  // unpack drug concentration data
  Rcpp::List data_list = args_params["data"];
  Rcpp::List data_drug = data_list["data_drug"];
  vector<int> ind = rcpp_to_vector_int(data_drug["individual"]);
  vector<int> time = rcpp_to_vector_int(data_drug["time"]);
  vector<double> drug_conc_raw = rcpp_to_vector_double(data_drug["drug_value.conc"]);
  
  // unpack data on individual group weights
  ind_weight = rcpp_to_vector_double(data_list["ind_weight"]);
  n_ind = ind_weight.size();
  n_time = drug_conc_raw.size() / n_ind;
  
  // unpack trial data
  Rcpp::List data_control = data_list["data_control"];
  control_n = rcpp_to_vector_int(data_control["n_patients"]);
  control_n_inf = rcpp_to_vector_int(data_control["n_infected"]);
  control_time0 = rcpp_to_vector_int(data_control["time"]);
  control_time1 = rcpp_to_vector_int(data_control["time.1"]);
  
  Rcpp::List data_treat = data_list["data_treat"];
  treat_n = rcpp_to_vector_int(data_treat["n_patients"]);
  treat_n_inf = rcpp_to_vector_int(data_treat["n_infected"]);
  treat_time0 = rcpp_to_vector_int(data_treat["time"]);
  treat_time1 = rcpp_to_vector_int(data_treat["time.1"]);
  
  // reformat drug concentration data into matrix
  drug_conc = vector<vector<double>>(n_ind, vector<double>(n_time));
  for (size_t i = 0; i < ind.size(); ++i) {
    drug_conc[ind[i] - 1][time[i]] = drug_conc_raw[i];
  }
}
