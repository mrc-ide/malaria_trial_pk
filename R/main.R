
#------------------------------------------------
#' @title Run drjacoby MCMC
#'
#' @description Run MCMC either with or without parallel tempering turned on.
#'   Minimum inputs include a data object, a data.frame of parameters, a
#'   log-likelihood function and a log-prior function. Produces an object of
#'   class \code{drjacoby_output}, which contains all MCMC output along with
#'   some diagnostics and a record of inputs.
#'   
#' @details Note that both \code{data} and \code{misc} are passed into
#'   log-likelihood/log-prior functions *by reference*. This means if you modify
#'   these objects inside the functions then any changes will persist.
#'
#' @param data a named list of numeric data values. When using C++ likelihood
#'   and/or prior these values are treated internally as doubles, so while
#'   integer and Boolean values can be used, keep in mind that these will be
#'   recast as doubles in the likelihood (i.e. \code{TRUE = 1.0}).
#' @param burnin the number of burn-in iterations. Automatic tuning of proposal
#'   standard deviations is only active during the burn-in period.
#' @param samples the number of sampling iterations.
#' @param rungs the number of temperature rungs used in the parallel tempering
#'   method. By default, \eqn{\beta} values are equally spaced between 0 and 1,
#'   i.e. \eqn{\beta[i]=}\code{(i-1)/(rungs-1)} for \code{i} in \code{1:rungs}.
#'   The likelihood for the \out{i<sup>th</sup>} heated chain is raised to the
#'   power \eqn{\beta[i]^\alpha}, meaning we can use the \eqn{\alpha} parameter
#'   to concentrate rungs towards the start or the end of the interval (see the
#'   \code{alpha} argument).
#' @param chains the number of independent replicates of the MCMC to run. If a
#'   \code{cluster} object is defined then these chains are run in parallel,
#'   otherwise they are run in serial.
#' @param beta_manual vector of manually defined \eqn{\beta} values used in the
#'   parallel tempering approach. If defined, this overrides the spacing defined
#'   by \code{rungs}. Note that even manually defined \eqn{\beta} values are
#'   raised to the power \eqn{\alpha} internally, hence you should set
#'   \code{alpha = 1} if you want to fix \eqn{\beta} values exactly.
#' @param alpha the likelihood for the \out{i<sup>th</sup>} heated chain is
#'   raised to the power \eqn{\beta[i]^\alpha}, meaning we can use the
#'   \eqn{\alpha} parameter to concentrate rungs towards the start or the end of
#'   the temperature scale.
#' @param target_acceptance Target acceptance rate. Should be between 0 and 1.
#'   Default of 0.44, set as optimum for unvariate proposal distributions.
#' @param cluster option to pass in a cluster environment, allowing chains to be
#'   run in parallel (see package "parallel").
#' @param coupling_on whether to implement Metropolis-coupling over temperature
#'   rungs. The option of deactivating coupling has been retained for general
#'   interest and debugging purposes only. If this parameter is \code{FALSE}
#'   then parallel tempering will have no impact on MCMC mixing.
#' @param pb_markdown whether to run progress bars in markdown mode, meaning
#'   they are only updated when they reach 100\% to avoid large amounts of output
#'   being printed to markdown files.
#' @param silent whether to suppress all console output.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats setNames var runif
#' @export

run_mcmc <- function(data,
                     burnin = 1e3,
                     samples = 1e4,
                     chains = 5,
                     beta_manual = 1,
                     target_acceptance = 0.44,
                     cluster = NULL,
                     coupling_on = TRUE,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  # ---------- check inputs ----------
  
  # check MCMC parameters
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  assert_vector_bounded(beta_manual)
  assert_increasing(beta_manual)
  rungs <- length(beta_manual)
  assert_eq(beta_manual[rungs], 1.0)
  assert_single_logical(coupling_on)
  assert_bounded(target_acceptance, 0, 1)
  
  # check misc parameters
  if (!is.null(cluster)) {
    assert_class(cluster, "cluster")
  }
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(data = data,
                      burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      coupling_on = coupling_on,
                      beta = beta_manual,
                      pb_markdown = pb_markdown,
                      silent = silent,
                      target_acceptance = target_acceptance)
  
  # functions to pass to C++
  args_functions <- list(update_progress = update_progress)
  
  # complete list of arguments
  args <- list(args_params = args_params,
               args_functions = args_functions)
  
  # create distinct argument sets over chains
  parallel_args <- replicate(chains, args, simplify = FALSE)
  for (i in 1:chains) {
    parallel_args[[i]]$args_params$chain <- i
  }
  
  
  # ---------- run MCMC ----------
  
  # split into parallel and serial implementations
  if (!is.null(cluster)) {
    
    # run in parallel
    parallel::clusterEvalQ(cluster, library(drjacoby))
    output_raw <- parallel::clusterApplyLB(cl = cluster, parallel_args, deploy_chain)
    
  } else {
    
    # run in serial
    output_raw <- lapply(parallel_args, deploy_chain)
  }
  
  # print total runtime
  chain_runtimes <- mapply(function(x) x$t_diff, output_raw)
  if (!silent) {
    message(sprintf("total MCMC run-time: %s seconds", signif(sum(chain_runtimes), 3)))
  }
  
  #return(output_raw)
  
  
  # ---------- process output ----------
  
  # make main output data.frame
  df_output <- mapply(function(i) {
    data.frame(chain = i,
               phase = c(rep("burnin", burnin), rep("sampling", samples)),
               iteration = 1:(burnin + samples),
               lambda = output_raw[[i]]$lambda,
               min_prob = output_raw[[i]]$min_prob,
               half_point = output_raw[[i]]$half_point,
               hill_power = output_raw[[i]]$hill_power,
               logprior = output_raw[[i]]$logprior,
               loglikelihood = output_raw[[i]]$loglike)
  }, seq_along(output_raw), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(chain = as.factor(chain))
  
  return(df_output)
  
  # append to output list
  output_processed <- list(output = df_output,
                           pt = df_pt)
  
  ## Diagnostics
  output_processed$diagnostics <- list()
  
  # run-times
  run_time <- data.frame(chain = chain_names,
                         seconds = chain_runtimes)
  output_processed$diagnostics$run_time <- run_time
  
  # Rhat (Gelman-Rubin diagnostic)
  if (chains > 1) {
    rhat_est <- c()
    for (p in seq_along(param_names)) {
      rhat_est[p] <- df_output %>%
        dplyr::filter(phase == "sampling") %>%
        dplyr::select(chain, param_names[p]) %>%
        gelman_rubin(chains = chains, samples = samples)
    }
    rhat_est[skip_param] <- NA
    names(rhat_est) <- param_names
    output_processed$diagnostics$rhat <- rhat_est
  }
  
  # ESS
  ess_est <- df_output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::select(param_names) %>%
    apply(2, coda::effectiveSize)
  ess_est[skip_param] <- NA
  output_processed$diagnostics$ess <- ess_est
  
  # Thermodynamic power
  output_processed$diagnostics$rung_details <- data.frame(rung = 1:rungs,
                                                          thermodynamic_power = beta_manual)
  
  # Metropolis-coupling
  # store acceptance rates between pairs of rungs (links)
  mc_accept <- NA
  if (rungs > 1) {
    
    # MC accept
    mc_accept <- expand.grid(link = seq_len(rungs - 1), chain = chain_names)
    mc_accept_burnin <- unlist(lapply(output_raw, function(x){x$mc_accept_burnin})) / burnin
    mc_accept_sampling <- unlist(lapply(output_raw, function(x){x$mc_accept_sampling})) / samples
    mc_accept <- rbind(cbind(mc_accept, phase = "burnin", value = mc_accept_burnin),
                       cbind(mc_accept, phase = "sampling", value = mc_accept_sampling))
    
  }
  output_processed$diagnostics$mc_accept <- mc_accept
  
  # DIC
  DIC <- df_pt %>%
    dplyr::filter(.data$phase == "sampling" & .data$rung == rungs) %>%
    dplyr::select(.data$loglikelihood) %>%
    dplyr::mutate(deviance = -2*.data$loglikelihood) %>%
    dplyr::summarise(DIC = mean(.data$deviance) + 0.5*var(.data$deviance)) %>%
    dplyr::pull(.data$DIC)
  output_processed$diagnostics$DIC_Gelman <- DIC
  
  # save output as custom class
  class(output_processed) <- "drjacoby_output"
  
  # return
  return(output_processed)
}

#------------------------------------------------
# deploy main_mcmc for this chain
#' @noRd
deploy_chain <- function(args) {
  
  # get parameters
  burnin <- args$args_params$burnin
  samples <- args$args_params$samples
  
  # make progress bars
  pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args$args_progress <- list(pb_burnin = pb_burnin,
                             pb_samples = pb_samples)
  
  # run C++ function
  ret <- main_cpp(args)
  
  # remove arguments
  rm(args)
  
  return(ret)
}

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

# Deal with user input cpp not being defined
globalVariables(c("create_xptr"))
