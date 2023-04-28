
#include "Particle.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s, double beta) {
  
  // pointer to system object
  this->s_ptr = &s;
  
  // thermodynamic power
  this->beta = beta;
  
  // parameters
  lambda = vector<double>(s_ptr->n_weeks, 0.02);
  min_prob = 0.5;
  half_point = 1.0;
  hill_power = 1.5;
  
  // random initial conditions
  for (size_t i = 0; i < lambda.size(); ++i) {
    lambda[i] = R::runif(0.0, 1.0);
  }
  min_prob = R::runif(0.0, 1.0);
  half_point = R::runif(0.0, 5.0);
  hill_power = R::runif(1.0, 5.0);
  
  // intermediate objects
  drug_pow = vector<vector<double>>(s_ptr->n_ind, vector<double>(s_ptr->n_time));
  drug_pow_prop = drug_pow;
  exp_rate = vector<vector<double>>(s_ptr->n_ind, vector<double>(s_ptr->treat_n.size()));
  
  // populate intermediate objects
  recalc_drug_pow(drug_pow, hill_power);
  
  // proposal parameters
  bw_lambda = 1.0;
  bw_min_prob = 1.0;
  bw_half_point = 1.0;
  bw_hill_power = 1.0;
  bw_stepsize = 1.0;
  
  // likelihoods and priors
  loglike_control = get_loglike_control(lambda);
  loglike_treat = get_loglike_treat(lambda, min_prob, half_point, hill_power, drug_pow);
  loglike = loglike_control + loglike_treat;
  logprior = get_logprior(lambda, min_prob, half_point, hill_power);
  
  // acceptance rates
  accept_count = 0;
}

//------------------------------------------------
void Particle::update(bool RM_on, int iteration) {
  update_lambda(RM_on, iteration);
  update_min_prob(RM_on, iteration);
  update_half_point(RM_on, iteration);
  update_hill_power(RM_on, iteration);
}

//------------------------------------------------
void Particle::update_lambda(bool RM_on, int iteration) {
  
  // copy lambda values
  vector<double> lambda_prop = lambda;
  
  // update all lambdas
  for (int i = 0; i < s_ptr->n_weeks; ++i) {
    
    // propose new value
    lambda_prop[i] = abs(R::rnorm(lambda[i], bw_lambda));
    
    // calculate likelihood and prior of proposed value
    double loglike_control_prop = get_loglike_control(lambda_prop);
    double loglike_treat_prop = get_loglike_treat(lambda_prop, min_prob, half_point, hill_power, drug_pow);
    double loglike_prop = loglike_control_prop + loglike_treat_prop;
    double logprior_prop = get_logprior(lambda_prop, min_prob, half_point, hill_power);
    
    // calculate Metropolis-Hastings ratio
    double MH;
    if (beta == 0.0){
      MH = (logprior_prop - logprior);
    } else {
      MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
    }
    
    // accept or reject move
    bool MH_accept = (log(R::runif(0,1)) < MH);
    
    // implement changes
    if (MH_accept) {
      
      // update parameters
      lambda[i] = lambda_prop[i];
      
      // update likelihoods
      loglike_control = loglike_control_prop;
      loglike_treat = loglike_treat_prop;
      loglike = loglike_prop;
      logprior = logprior_prop;
      
      // Robbins-Monro positive update  (on the log scale)
      if (RM_on) {
        bw_lambda = exp(log(bw_lambda) + bw_stepsize*(1 - s_ptr->target_acceptance) / sqrt(iteration));
      }
      
      // add to acceptance rate count
      accept_count++;
      
    } else {
      
      // update parameters
      lambda_prop[i] = lambda[i];
      
      // Robbins-Monro negative update (on the log scale)
      if (RM_on) {
        bw_lambda = exp(log(bw_lambda) - bw_stepsize*s_ptr->target_acceptance / sqrt(iteration));
      }
      
    } // end MH step
    
  }
}

//------------------------------------------------
void Particle::update_min_prob(bool RM_on, int iteration) {
  
  // propose new value
  double min_prob_prop = rnorm1_interval(min_prob, bw_min_prob);
  
  // calculate likelihood and prior of proposed value
  double loglike_treat_prop = get_loglike_treat(lambda, min_prob_prop, half_point, hill_power, drug_pow);
  double loglike_prop = loglike_control + loglike_treat_prop;
  double logprior_prop = get_logprior(lambda, min_prob_prop, half_point, hill_power);
  
  // calculate Metropolis-Hastings ratio
  double MH;
  if (beta == 0.0){
    MH = (logprior_prop - logprior);
  } else {
    MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  }
  
  // accept or reject move
  bool MH_accept = (log(R::runif(0,1)) < MH);
  
  // implement changes
  if (MH_accept) {
    
    // update parameters
    min_prob = min_prob_prop;
    
    // update likelihoods
    loglike_treat = loglike_treat_prop;
    loglike = loglike_prop;
    logprior = logprior_prop;
    
    // Robbins-Monro positive update  (on the log scale)
    if (RM_on) {
      bw_min_prob = exp(log(bw_min_prob) + bw_stepsize*(1 - s_ptr->target_acceptance) / sqrt(iteration));
    }
    
    // add to acceptance rate count
    accept_count++;
    
  } else {
    
    // Robbins-Monro negative update (on the log scale)
    if (RM_on) {
      bw_min_prob = exp(log(bw_min_prob) - bw_stepsize*s_ptr->target_acceptance / sqrt(iteration));
    }
    
  } // end MH step
  
}

//------------------------------------------------
void Particle::update_half_point(bool RM_on, int iteration) {
  
  // propose new value
  double half_point_prop = abs(R::rnorm(half_point, bw_half_point));
  
  // calculate likelihood and prior of proposed value
  double loglike_treat_prop = get_loglike_treat(lambda, min_prob, half_point_prop, hill_power, drug_pow);
  double loglike_prop = loglike_control + loglike_treat_prop;
  double logprior_prop = get_logprior(lambda, min_prob, half_point_prop, hill_power);
  
  // calculate Metropolis-Hastings ratio
  double MH;
  if (beta == 0.0){
    MH = (logprior_prop - logprior);
  } else {
    MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  }
  
  // accept or reject move
  bool MH_accept = (log(R::runif(0,1)) < MH);
  
  // implement changes
  if (MH_accept) {
    
    // update parameters
    half_point = half_point_prop;
    
    // update likelihoods
    loglike_treat = loglike_treat_prop;
    loglike = loglike_prop;
    logprior = logprior_prop;
    
    // Robbins-Monro positive update  (on the log scale)
    if (RM_on) {
      bw_half_point = exp(log(bw_half_point) + bw_stepsize*(1 - s_ptr->target_acceptance) / sqrt(iteration));
    }
    
    // add to acceptance rate count
    accept_count++;
    
  } else {
    
    // Robbins-Monro negative update (on the log scale)
    if (RM_on) {
      bw_half_point = exp(log(bw_half_point) - bw_stepsize*s_ptr->target_acceptance / sqrt(iteration));
    }
    
  } // end MH step
  
}

//------------------------------------------------
void Particle::update_hill_power(bool RM_on, int iteration) {
  
  // propose new value
  double hill_power_prop = abs(R::rnorm(hill_power, bw_hill_power));
  //double hill_power_prop = rnorm1_interval(hill_power, bw_hill_power, 1.0, 20.0);
  
  // recalculate intermediate values
  recalc_drug_pow(drug_pow_prop, hill_power_prop);
  
  // calculate likelihood and prior of proposed value
  double loglike_treat_prop = get_loglike_treat(lambda, min_prob, half_point, hill_power_prop, drug_pow_prop);
  double loglike_prop = loglike_control + loglike_treat_prop;
  double logprior_prop = get_logprior(lambda, min_prob, half_point, hill_power_prop);
  
  // calculate Metropolis-Hastings ratio
  double MH;
  if (beta == 0.0){
    MH = (logprior_prop - logprior);
  } else {
    MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior);
  }
  
  // accept or reject move
  bool MH_accept = (log(R::runif(0,1)) < MH);
  
  // implement changes
  if (MH_accept) {
    
    // update parameters
    hill_power = hill_power_prop;
    
    // update intermediate values
    drug_pow = drug_pow_prop;
    
    // update likelihoods
    loglike_treat = loglike_treat_prop;
    loglike = loglike_prop;
    logprior = logprior_prop;
    
    // Robbins-Monro positive update  (on the log scale)
    if (RM_on) {
      bw_hill_power = exp(log(bw_hill_power) + bw_stepsize*(1 - s_ptr->target_acceptance) / sqrt(iteration));
    }
    
    // add to acceptance rate count
    accept_count++;
    
  } else {
    
    // Robbins-Monro negative update (on the log scale)
    if (RM_on) {
      bw_hill_power = exp(log(bw_hill_power) - bw_stepsize*s_ptr->target_acceptance / sqrt(iteration));
    }
    
  } // end MH step
  
}

//------------------------------------------------
// recalculate matrix based on power parameter of hill function
void Particle::recalc_drug_pow(vector<vector<double>> &mat, double k) {
  //return;
  
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    for (int j = 0; j < s_ptr->n_time; ++j) {
      mat[i][j] = pow(s_ptr->drug_conc[i][j], k);
    }
  }
}

//------------------------------------------------
double Particle::get_loglike_control(vector<double> &lambda_) {
  //return 0.0;
  
  double ret = 0.0;
  for (size_t i = 0; i < s_ptr->control_n.size(); ++i) {
    double r = 0.0;
    for (size_t j = 0; j < s_ptr->control_lambda_index[i].size(); ++j) {
      int lambda_index = s_ptr->control_lambda_index[i][j] - 1;
      r += lambda_[lambda_index] / 24 * s_ptr->control_lambda_weight[i][j];
    }
    double p_sus = exp(-r);
    ret += R::dbinom(s_ptr->control_n_inf[i], s_ptr->control_n[i], 1.0 - p_sus, true);
  }
  
  return ret;
}

//------------------------------------------------
double Particle::get_loglike_treat(vector<double> &lambda_, double min_prob_, double half_point_, 
                                   double hill_power_, vector<vector<double>> &drug_pow_) {
  //return 0.0;
  
  double h_raised = pow(half_point_, hill_power_);
  
  // calculate the rate of the exponential function for each trial window
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    fill(exp_rate[i].begin(), exp_rate[i].end(), 0.0);
    int window = 0;
    int week_index = 0;
    for (int j = 0; j < s_ptr->n_time; ++j) {
      
      // move week window forward if needed
      if (j == (week_index + 1)*7*24) {
        week_index++;
        if (week_index >= lambda_.size()) {
          week_index = lambda_.size() - 1;
        }
      }
      
      // move window forward if needed, or break if reached end
      if (j == s_ptr->treat_time1[window]) {
        window++;
        if (window == s_ptr->treat_n.size()) {
          break;
        }
      }
      
      double hill = min_prob_ + (1.0 - min_prob_) * h_raised / (h_raised + drug_pow_[i][j]);
      exp_rate[i][window] += lambda_[week_index] / 24 * hill;
    }
  }
  
  // calculate likelihood over trial data
  double ret = 0.0;
  for (size_t t = 0; t < s_ptr->treat_n.size(); ++t) {
    double p_sus = 0.0;
    for (int i = 0; i < s_ptr->n_ind; ++i) {
      p_sus += s_ptr->ind_weight[i] * exp(-exp_rate[i][t]);
    }
    ret += R::dbinom(s_ptr->treat_n_inf[t], s_ptr->treat_n[t], 1.0 - p_sus, true);
  }
  
  return ret;
}

//------------------------------------------------
double Particle::get_logprior(vector<double> &lambda_, double min_prob_,
                              double half_point_, double hill_power_) {
  
  return 0.0;
  
  // lognormal prior on lambdas
  double ret = 0.0;
  double lambda_mean = 0.01;
  double lambda_sd = 0.05;
  double lambda_sigma2 = log(pow(lambda_sd / lambda_mean, 2.0) + 1.0);
  double lambda_sigma = pow(lambda_sigma2, 0.5);
  double lambda_mu = log(lambda_mean) - lambda_sigma2 / 2.0;
  for (size_t i = 0; i < lambda_.size(); ++i) {
    ret += R::dlnorm(lambda_[i], lambda_mu, lambda_sigma, true);
  }
  
  // beta prior on min_prob
  ret += R::dbeta(min_prob_, 1.0, 19.0, true);
  
  // lognormal prior on half_point
  //double half_point_mean = 1.5;
  //double half_point_sd = 1.0;
  //double half_point_sigma2 = log(pow(half_point_sd / half_point_mean, 2.0) + 1.0);
  //double half_point_mu = log(half_point_mean) - half_point_sigma2 / 2.0;
  //ret += R::dlnorm(half_point_, half_point_mu, pow(half_point_sigma2, 0.5), true);
  
  // half-normal prior on half_point
  ret += R::dnorm4(half_point_, 0.0, 5.0, true);
  
  // half-normal prior on hill_power
  ret += R::dnorm4(hill_power_, 0.0, 5.0, true);
  
  return ret;
}

//------------------------------------------------
double Particle::get_loglike_fromparams(Rcpp::List params) {
  
  // import parameter values
  for (size_t i = 0; i < lambda.size(); ++i) {
    lambda[i] = params[i];
  }
  min_prob = params["min_prob"];
  half_point = params["half_point"];
  hill_power = params["hill_power"];
  
  // populate intermediate objects
  recalc_drug_pow(drug_pow, hill_power);
  
  // calculate likelihood
  loglike_control = get_loglike_control(lambda);
  loglike_treat = get_loglike_treat(lambda, min_prob, half_point, hill_power, drug_pow);
  loglike = loglike_control + loglike_treat;
  
  return loglike;
}
