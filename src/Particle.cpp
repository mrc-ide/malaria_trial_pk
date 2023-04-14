
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
  lambda = 1.0;
  min_prob = 0.5;
  half_point = 1.0;
  hill_power = 1.5;
  
  // intermediate objects
  drug_pow = vector<vector<double>>(s_ptr->n_ind, vector<double>(s_ptr->n_time));
  drug_pow_prop = drug_pow;
  sum_hill_unscaled = vector<vector<double>>(s_ptr->n_ind, vector<double>(s_ptr->treat_n.size()));
  sum_hill_unscaled_prop = sum_hill_unscaled;
  
  // populate intermediate objects
  recalc_drug_pow(drug_pow, hill_power);
  recalc_sum_hill_unscaled(sum_hill_unscaled, half_point, drug_pow);
  
  // proposal parameters
  bw_lambda = 1.0;
  bw_min_prob = 1.0;
  bw_half_point = 1.0;
  bw_hill_power = 1.0;
  bw_stepsize = 1.0;
  
  // likelihoods and priors
  loglike_control = get_loglike_control(lambda);
  loglike_treat = get_loglike_treat(lambda, min_prob, sum_hill_unscaled);
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
  
  // propose new value
  double lambda_prop = abs(R::rnorm(lambda, bw_lambda));
  
  // calculate likelihood and prior of proposed value
  double loglike_control_prop = get_loglike_control(lambda_prop);
  double loglike_treat_prop = get_loglike_treat(lambda_prop, min_prob, sum_hill_unscaled);
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
    lambda = lambda_prop;
    
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
    
    // Robbins-Monro negative update (on the log scale)
    if (RM_on) {
      bw_lambda = exp(log(bw_lambda) - bw_stepsize*s_ptr->target_acceptance / sqrt(iteration));
    }
    
  } // end MH step
}

//------------------------------------------------
void Particle::update_min_prob(bool RM_on, int iteration) {
  
  // propose new value
  double min_prob_prop = rnorm1_interval(min_prob, bw_min_prob);
  
  // calculate likelihood and prior of proposed value
  double loglike_treat_prop = get_loglike_treat(lambda, min_prob_prop, sum_hill_unscaled);
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
  
  // recalculate intermediate values
  recalc_sum_hill_unscaled(sum_hill_unscaled_prop, half_point_prop, drug_pow);
  
  // calculate likelihood and prior of proposed value
  double loglike_treat_prop = get_loglike_treat(lambda, min_prob, sum_hill_unscaled_prop);
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
    
    // update intermediate values
    sum_hill_unscaled = sum_hill_unscaled_prop;
    
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
  recalc_sum_hill_unscaled(sum_hill_unscaled_prop, half_point, drug_pow_prop);
  
  // calculate likelihood and prior of proposed value
  double loglike_treat_prop = get_loglike_treat(lambda, min_prob, sum_hill_unscaled_prop);
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
    sum_hill_unscaled = sum_hill_unscaled_prop;
    
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
// recalculate matrix based on half_point parameter of hill function
void Particle::recalc_sum_hill_unscaled(vector<vector<double>> &mat, double h,
                                        vector<vector<double>> &drug_pow_) {
  //return;
  
  double h_raised = pow(h, hill_power);
  
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    fill(mat[i].begin(), mat[i].end(), 0.0);
    int window = 0;
    for (int j = 0; j < s_ptr->n_time; ++j) {
      
      // move window forward if needed, or break if reached end
      if (j == s_ptr->treat_time1[window]) {
        window++;
        if (window == s_ptr->treat_n.size()) {
          break;
        }
      }
      
      mat[i][window] += h_raised / (h_raised + drug_pow_[i][j]);
    }
  }
}

//------------------------------------------------
double Particle::get_loglike_control(double lambda_) {
  //return 0.0;
  
  double lambda_hourly = lambda_ / 24;
  
  double ret = 0.0;
  for (size_t i = 0; i < s_ptr->control_n.size(); ++i) {
    double t_diff = s_ptr->control_time1[i] - s_ptr->control_time0[i];
    double p_sus = exp(-lambda_hourly*t_diff);
    ret += R::dbinom(s_ptr->control_n_inf[i], s_ptr->control_n[i], 1.0 - p_sus, true);
  }
  
  return ret;
}

//------------------------------------------------
double Particle::get_loglike_treat(double lambda_, double min_prob_,
                                   vector<vector<double>> &sum_hill_unscaled_) {
  //return 0.0;
  
  double lambda_hourly = lambda_ / 24;
  
  double ret = 0.0;
  for (size_t t = 0; t < s_ptr->treat_n.size(); ++t) {
    double t_diff = s_ptr->treat_time1[t] - s_ptr->treat_time0[t];
    double p_sus = 0.0;
    for (int i = 0; i < s_ptr->n_ind; ++i) {
      p_sus += s_ptr->ind_weight[i] * exp(-lambda_hourly * (t_diff*min_prob_ + (1.0 - min_prob_)*sum_hill_unscaled_[i][t]));
    }
    ret += R::dbinom(s_ptr->treat_n_inf[t], s_ptr->treat_n[t], 1.0 - p_sus, true);
  }
  
  return ret;
}

//------------------------------------------------
double Particle::get_logprior(double lambda_, double min_prob_,
                              double half_point_, double hill_power_) {
  
  // lognormal prior on lambda
  double ret = 0.0;
  double lambda_mean = 0.01;
  double lambda_sd = 0.05;
  double lambda_sigma2 = log(pow(lambda_sd / lambda_mean, 2.0) + 1.0);
  double lambda_mu = log(lambda_mean) - lambda_sigma2 / 2.0;
  ret += R::dlnorm(lambda_, lambda_mu, pow(lambda_sigma2, 0.5), true);
  
  // beta prior on min_prob
  ret += R::dbeta(min_prob_, 1.0, 19.0, true);
  
  // lognormal prior on half-point
  double half_point_mean = 1.5;
  double half_point_sd = 1.0;
  double half_point_sigma2 = log(pow(half_point_sd / half_point_mean, 2.0) + 1.0);
  double half_point_mu = log(half_point_mean) - half_point_sigma2 / 2.0;
  ret += R::dlnorm(half_point_, half_point_mu, pow(half_point_sigma2, 0.5), true);
  
  // half-normal prior on hill_power
  ret += R::dnorm4(hill_power_, 0.0, 2.0, true);
  
  return ret;
}
