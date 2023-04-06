
#include "System.h"
#include "misc.h"

using namespace std;

void System::load(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
  // misc
  misc = args_params["misc"];
  
  // model parameters
  theta_vector = args_params["theta_vector"];
  theta_min = rcpp_to_vector_double(args_params["theta_min"]);
  theta_max = rcpp_to_vector_double(args_params["theta_max"]);
  block = rcpp_to_matrix_int(args_params["block"]);
  n_block = rcpp_to_int(args_params["n_block"]);
  trans_type = rcpp_to_vector_int(args_params["trans_type"]);
  skip_param = rcpp_to_vector_bool(args_params["skip_param"]);
  d = int(theta_min.size());
  target_acceptance = rcpp_to_double(args_params["target_acceptance"]);
  
  // MCMC parameters
  burnin = rcpp_to_int(args_params["burnin"]);
  samples = rcpp_to_int(args_params["samples"]);
  rungs = rcpp_to_int(args_params["rungs"]);
  coupling_on = rcpp_to_bool(args_params["coupling_on"]);
  beta_raised = rcpp_to_vector_double(args_params["beta_raised"]);
  chain = rcpp_to_int(args_params["chain"]);
  
  // misc parameters
  save_hot_draws = rcpp_to_bool(args_params["save_hot_draws"]);
  pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  silent = rcpp_to_bool(args_params["silent"]);
  
  
  // unpack raw data
  if (!silent) {
    print("Loading data");
  }
  Rcpp::List data_list = args_params["data"];
  vector<int> ind = rcpp_to_vector_int(data_list["ind"]);
  vector<int> time = rcpp_to_vector_int(data_list["time"]);
  vector<double> drug_conc_raw = rcpp_to_vector_double(data_list["drug_conc"]);
  ind_weight = rcpp_to_vector_double(data_list["ind_weight"]);
  n_ind = ind_weight.size();
  n_time = drug_conc_raw.size() / n_ind;
  
  // reformat into matrix
  drug_conc = vector<vector<double>>(n_ind, vector<double>(n_time));
  for (int i = 0; i < n_ind; ++i) {
    drug_conc[ind[i] - 1][time[i]] = drug_conc_raw[i];
  }
  
}
