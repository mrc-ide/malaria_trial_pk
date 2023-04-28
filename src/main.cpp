
#include "main.h"
#include "misc.h"
#include "System.h"
#include <chrono>

using namespace std;

//------------------------------------------------
// main MCMC function
Rcpp::List main_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = chrono::high_resolution_clock::now();
  
  // create sytem object and load args
  System s;
  s.load_data(args["args_data"]);
  s.load_params(args["args_params"]);
  
  // extract R utility functions that will be called from within MCMC
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // extract progress bar objects
  Rcpp::List args_progress = args["args_progress"];
  
  // local copies of some parameters for convenience
  int rungs = s.rungs;
  int cold_chain = rungs - 1;
  
  // initialise vector of particles
  vector<Particle> particle_vec(rungs);
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].init(s, s.beta[r]);
  }
  
  // objects for storing loglikelihood parameters
  vector<double> loglike(s.burnin + s.samples);
  vector<double> logprior(s.burnin + s.samples);
  vector<vector<double>> lambda(s.burnin + s.samples);
  vector<double> min_prob(s.burnin + s.samples);
  vector<double> half_point(s.burnin + s.samples);
  vector<double> hill_power(s.burnin + s.samples);
  
  // specify stored values at first iteration
  loglike[0] = particle_vec[cold_chain].loglike;
  logprior[0] = particle_vec[cold_chain].logprior;
  lambda[0] = particle_vec[cold_chain].lambda;
  min_prob[0] = particle_vec[cold_chain].min_prob;
  half_point[0] = particle_vec[cold_chain].half_point;
  hill_power[0] = particle_vec[cold_chain].hill_power;
  
  // store Metropolis coupling acceptance rates
  vector<int> mc_accept_burnin(rungs - 1);
  vector<int> mc_accept_sampling(rungs - 1);
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("MCMC chain", s.chain);
    print("burn-in");
  }
  
  // loop through burn-in iterations
  for (int rep = 1; rep < s.burnin; ++rep) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // loop through rungs and update particles
    for (int r = 0; r < rungs; ++r) {
      particle_vec[r].update(true, rep);
    }
    
    // store results
    loglike[rep] = particle_vec[cold_chain].loglike;
    logprior[rep] = particle_vec[cold_chain].logprior;
    lambda[rep] = particle_vec[cold_chain].lambda;
    min_prob[rep] = particle_vec[cold_chain].min_prob;
    half_point[rep] = particle_vec[cold_chain].half_point;
    hill_power[rep] = particle_vec[cold_chain].hill_power;
    
    // perform Metropolis coupling
    if (s.coupling_on) {
      coupling(particle_vec, mc_accept_burnin);
    }
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.burnin) / 100));
      if ((remainder == 0 && !s.pb_markdown) || ((rep+1) == s.burnin)) {
        update_progress(args_progress, "pb_burnin", rep+1, s.burnin, false);
        if ((rep + 1) == s.burnin) {
          print("");
        }
      }
    }
    
  }  // end burn-in MCMC loop
  
  // print phase diagnostics
  if (!s.silent) {
    //double accept_rate = particle_vec[cold_chain].accept_count / double(s.burnin*d);
    //Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000) / 10.0 << "%\n";
  }
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("sampling phase");
  }
  
  // reset acceptance count of all rungs
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].accept_count = 0;
  }
  
  // loop through sampling iterations
  for (int rep = 0; rep < s.samples; ++rep) {
    
    // allow user to exit on escape
    Rcpp::checkUserInterrupt();
    
    // loop through rungs and update particles
    for (int r = 0; r < rungs; ++r) {
      particle_vec[r].update(false, rep);
    }
    
    // store results
    loglike[s.burnin + rep] = particle_vec[cold_chain].loglike;
    logprior[s.burnin + rep] = particle_vec[cold_chain].logprior;
    lambda[s.burnin + rep] = particle_vec[cold_chain].lambda;
    min_prob[s.burnin + rep] = particle_vec[cold_chain].min_prob;
    half_point[s.burnin + rep] = particle_vec[cold_chain].half_point;
    hill_power[s.burnin + rep] = particle_vec[cold_chain].hill_power;
    
    // perform Metropolis coupling
    if (s.coupling_on) {
      coupling(particle_vec, mc_accept_sampling);
    }
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.samples) / 100));
      if ((remainder == 0 && !s.pb_markdown) || ((rep + 1) == s.samples)) {
        update_progress(args_progress, "pb_samples", rep + 1, s.samples, false);
        if ((rep + 1) == s.samples) {
          print("");
        }
      }
    }
    
  }  // end sampling MCMC loop
  
  // print final diagnostics
  if (!s.silent) {
    //double accept_rate = particle_vec[cold_chain].accept_count / double(s.samples*d);
    //Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000) / 10.0 << "%\n";
  }
  
  
  // ---------- return ----------
  
  // end timer
  double t_diff = chrono_timer(t0, "chain completed in ", !s.silent);
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike") = loglike,
                                      Rcpp::Named("logprior") = logprior,
                                      Rcpp::Named("lambda") = lambda,
                                      Rcpp::Named("min_prob") = min_prob,
                                      Rcpp::Named("half_point") = half_point,
                                      Rcpp::Named("hill_power") = hill_power,
                                      Rcpp::Named("mc_accept_burnin") = mc_accept_burnin,
                                      Rcpp::Named("mc_accept_sampling") = mc_accept_sampling,
                                      Rcpp::Named("t_diff") = t_diff);
  return ret;
}

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(vector<Particle> &particle_vec, vector<int> &mc_accept) {
  
  // get number of rungs
  int rungs = int(particle_vec.size());
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i = 0; i < (rungs-1); ++i) {
    
    // define rungs of interest
    int rung1 = i;
    int rung2 = i + 1;
    
    // get log-likelihoods and beta values of two chains in the comparison
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta1 = particle_vec[rung1].beta;
    double beta2 = particle_vec[rung2].beta;
    
    // calculate acceptance ratio (still in log space)
    double acceptance;
    if (beta1 == 0.0){
      acceptance = (loglike1*beta2) - (loglike2*beta2);
    } else {
      acceptance = (loglike2*beta1 + loglike1*beta2) - (loglike1*beta1 + loglike2*beta2);
    }
    
    // accept or reject move
    bool accept_move = (log(R::runif(0,1)) < acceptance);
    
    // implement swap
    if (accept_move) {
      
      // swap parameter values
      vector<double> tmp_vec = particle_vec[rung1].lambda;
      particle_vec[rung1].lambda = particle_vec[rung2].lambda;
      particle_vec[rung2].lambda = tmp_vec;
      
      double tmp = particle_vec[rung1].min_prob;
      particle_vec[rung1].min_prob = particle_vec[rung2].min_prob;
      particle_vec[rung2].min_prob = tmp;
      
      tmp = particle_vec[rung1].half_point;
      particle_vec[rung1].half_point = particle_vec[rung2].half_point;
      particle_vec[rung2].half_point = tmp;
      
      tmp = particle_vec[rung1].hill_power;
      particle_vec[rung1].hill_power = particle_vec[rung2].hill_power;
      particle_vec[rung2].hill_power = tmp;
      
      // swap loglikelihoods and logpriors
      tmp = particle_vec[rung1].loglike;
      particle_vec[rung1].loglike = particle_vec[rung2].loglike;
      particle_vec[rung2].loglike = tmp;
      
      tmp = particle_vec[rung1].logprior;
      particle_vec[rung1].logprior = particle_vec[rung2].logprior;
      particle_vec[rung2].logprior = tmp;
      
      // update acceptance rates
      mc_accept[i]++;
    }
    
  }  // end loop over rungs
}

//------------------------------------------------
// return loglikelihood
double get_loglike_cpp(Rcpp::List args) {
  
  // create sytem object and load data only
  System s;
  s.load_data(args["args_data"]);
  
  // create single particle and initialise
  Particle p;
  p.init(s, 1.0);
  
  // get likelihood
  double ret = p.get_loglike_fromparams(args["params"]);
  return ret;
}
