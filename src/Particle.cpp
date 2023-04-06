
#include "Particle.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s, double beta_raised) {
  
  // pointer to system object
  this->s_ptr = &s;
  
  // local copies of some parameters for convenience
  d = s_ptr->d;
  
  // thermodynamic power
  this->beta_raised = beta_raised;
  
  // theta is the parameter vector in natural space
  theta = Rcpp::clone(s_ptr->theta_vector);
  theta_prop = Rcpp::clone(theta);
  
  // phi is a vector of transformed parameters
  phi = Rcpp::NumericVector(d);
  theta_to_phi();
  phi_prop = Rcpp::NumericVector(d);
  
  // proposal parameters
  bw = vector<double>(d, 1.0);
  bw_index = vector<int>(d, 1);
  bw_stepsize = 1.0;
  
  // likelihoods and priors
  loglike_block = vector<double>(s_ptr->n_block);
  loglike_prop_block = vector<double>(s_ptr->n_block);
  loglike = 0;
  loglike_prop = 0;
  logprior = 0;
  logprior_prop = 0;
  
  // acceptance rates
  accept_count = 0;
  
}

//------------------------------------------------
// propose new value of phi[i] from univariate normal distribution
void Particle::propose_phi(int i) {
  phi_prop[i] = R::rnorm(phi[i], bw[i]);
}

//------------------------------------------------
// transform phi_prop to theta_prop. See main.R for a key to transformation
// types
void Particle::phi_prop_to_theta_prop(int i) {
  
  switch(s_ptr->trans_type[i]) {
  case 0:
    theta_prop[i] = phi_prop[i];
    break;
  case 1:
    theta_prop[i] = s_ptr->theta_max[i] - exp(phi_prop[i]);
    break;
  case 2:
    theta_prop[i] = exp(phi_prop[i]) + s_ptr->theta_min[i];
    break;
  case 3:
    theta_prop[i] = (s_ptr->theta_max[i]*exp(phi_prop[i]) + s_ptr->theta_min[i]) / (1 + exp(phi_prop[i]));
    break;
  default:
    Rcpp::stop("trans_type invalid");
  }
  
}

//------------------------------------------------
// transform theta to phi. See main.R for a key to transformation types
void Particle::theta_to_phi() {
  
  for (int i = 0; i < d; ++i) {
    switch(s_ptr->trans_type[i]) {
    case 0:
      phi[i] = theta[i];
      break;
    case 1:
      phi[i] = log(s_ptr->theta_max[i] - theta[i]);
      break;
    case 2:
      phi[i] = log(theta[i] - s_ptr->theta_min[i]);
      break;
    case 3:
      phi[i] = log(theta[i] - s_ptr->theta_min[i]) - log(s_ptr->theta_max[i] - theta[i]);
      break;
    default:
      Rcpp::stop("trans_type invalid");
    }
  }
  
}

//------------------------------------------------
// get adjustment factor to account for reparameterisation
double Particle::get_adjustment(int i) {
  
  double ret = 0;
  switch(s_ptr->trans_type[i]) {
  case 0:
    // (no adjustment needed)
    break;
  case 1:
    ret = log(s_ptr->theta_max[i] - theta_prop[i]) - log(s_ptr->theta_max[i] - theta[i]);
    break;
  case 2:
    ret = log(theta_prop[i] - s_ptr->theta_min[i]) - log(theta[i] - s_ptr->theta_min[i]);
    break;
  case 3:
    ret = log(s_ptr->theta_max[i] - theta_prop[i]) + log(theta_prop[i] - s_ptr->theta_min[i]) - log(s_ptr->theta_max[i] - theta[i]) - log(theta[i] - s_ptr->theta_min[i]);
    break;
  default:
    Rcpp::stop("trans_type invalid");
  }
  return ret;
  
}

//------------------------------------------------
// fixed log-likelihood
double Particle::get_loglike(Rcpp::NumericVector params) {
  
  // extract parameters
  double FOI = params["FOI"];
  double hourly_FOI = FOI / 24.0;
  double min_prob = params["min_prob"];
  double half_point = params["half_point"];
  double hill_power = params["hill_power"];
  
  // dummy example of how we can loop over the matrix of drug concentrations
  double ret = 0.0;
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    for (int j = 0; j < s_ptr->n_time; ++j) {
      double prob_susceptible = min_prob + (1.0 - min_prob) / (1.0 + pow(s_ptr->drug_conc[i][j] / half_point, hill_power));
      double prob_infection = 1.0 - exp(- hourly_FOI * prob_susceptible);
    }
  }
  
  return ret;
}

//------------------------------------------------
// fixed log-prior
double Particle::get_logprior(Rcpp::NumericVector params) {
  
  // extract parameters
  double FOI = params["FOI"];
  double min_prob = params["min_prob"];
  double half_point = params["half_point"];
  double hill_power = params["hill_power"];
  
  // lognormal prior on FOI
  double ret = 0.0;
  double FOI_mean = 0.01;
  double FOI_sd = 0.05;
  double FOI_sigma2 = log(pow(FOI_sd / FOI_mean, 2.0) + 1.0);
  double FOI_mu = log(FOI_mean) - FOI_sigma2 / 2.0;
  ret += R::dlnorm(FOI, FOI_mu, pow(FOI_sigma2, 0.5), true);
  
  // beta prior on min_prob
  ret += R::dbeta(min_prob, 1.0, 19.0, true);
  
  // lognormal prior on half-point
  double half_point_mean = 1.5;
  double half_point_sd = 1.0;
  double half_point_sigma2 = log(pow(half_point_sd / half_point_mean, 2.0) + 1.0);
  double half_point_mu = log(half_point_mean) - half_point_sigma2 / 2.0;
  ret += R::dlnorm(half_point, half_point_mu, pow(half_point_sigma2, 0.5), true);
  
  // half-normal prior on hill_power
  ret += R::dnorm4(hill_power, 0.0, 2.0, true);
  
  return ret;
}

