
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
#' @param chains the number of independent replicates of the MCMC to run.
#' @param beta_manual vector of manually defined \eqn{\beta} values used in the
#'   parallel tempering approach. If defined, this overrides the spacing defined
#'   by \code{rungs}. Note that even manually defined \eqn{\beta} values are
#'   raised to the power \eqn{\alpha} internally, hence you should set
#'   \code{alpha = 1} if you want to fix \eqn{\beta} values exactly.
#' @param target_acceptance Target acceptance rate. Should be between 0 and 1.
#'   Default of 0.44, set as optimum for unvariate proposal distributions.
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
                     coupling_on = TRUE,
                     pb_markdown = FALSE,
                     silent = FALSE) {
  
  # avoid no visible binding note
  chain <- iteration <- loglikelihood <- logprior <- phase <- NULL
  
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
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  # ---------- pre-processing ----------
  
  # get number of weeks in control data and define weekly breaks at hourly
  # resolution
  data_control <- data$data_control
  n_weeks <- ceiling(max(data_control$time.1) / 24 / 7)
  week_breaks <- (0:n_weeks)*7*24
  
  # populate two lists. control_lambda_index tells us which lambda parameter
  # values are needed for each window, and control_lambda_weight tells us the
  # corresponding weighting of these parameters. This is because some
  # observation windows can span multiple weeks, meaning a weighted sum of
  # lambda values is required
  control_lambda_index <- list()
  control_lambda_weight <- list()
  for (i in 1:nrow(data_control)) {
    window_days <- data_control$time[i]:(data_control$time.1[i] - 1)
    window_match <- as.numeric(cut(window_days, week_breaks, right = FALSE))
    control_lambda_index[[i]] <- unique(window_match)
    control_lambda_weight[[i]] <- as.vector(table(window_match))
  }
  
  # the eir_adjustment vector is needed in full for the treatment arm, but for
  # the control arm we just need the unique values in this vector and their
  # corresponding weights
  eir_unique <- unique(eir_adjustment)
  eir_weight <- mapply(sum, split(data$ind_weight, f = match(eir_adjustment, eir_unique)))
  
  # ---------- define argument lists ----------
  
  # data to pass to C++
  args_data <- c(data,
                 list(n_weeks = n_weeks,
                      control_lambda_index = control_lambda_index,
                      control_lambda_weight = control_lambda_weight,
                      eir_unique = eir_unique,
                      eir_weight = eir_weight))
  
  # parameters to pass to C++
  args_params <- list(burnin = burnin,
                      samples = samples,
                      rungs = rungs,
                      coupling_on = coupling_on,
                      beta = beta_manual,
                      pb_markdown = pb_markdown,
                      silent = silent,
                      target_acceptance = target_acceptance)
  
  # functions to pass to C++
  args_functions <- list(update_progress = update_progress)
  
  
  # ---------- run MCMC ----------
  
  # get output in list over chains
  output_raw <- list()
  for (i in 1:chains) {
    
    # index this chain
    args_params$chain <- i
    
    # make progress bars
    pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
    pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
    args_progress <- list(pb_burnin = pb_burnin,
                          pb_samples = pb_samples)
    
    # complete list of arguments
    args <- list(args_data = args_data,
                 args_params = args_params,
                 args_functions = args_functions,
                 args_progress = args_progress)
    
    # run C++ function
    output_raw[[i]] <- main_cpp(args)
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
    
    # process lambda values
    lambda_mat <- matrix(unlist(output_raw[[i]]$lambda), ncol = n_weeks, byrow = TRUE)
    colnames(lambda_mat) <- sprintf("lambda_%s", 1:n_weeks)
    
    # combine with other output
    data.frame(chain = i,
               phase = c(rep("burnin", burnin), rep("sampling", samples)),
               iteration = 1:(burnin + samples)) %>%
      bind_cols(lambda_mat) %>%
      bind_cols(data.frame(min_prob = output_raw[[i]]$min_prob,
                           half_point = output_raw[[i]]$half_point,
                           hill_power = output_raw[[i]]$hill_power,
                           logprior = output_raw[[i]]$logprior,
                           loglikelihood = output_raw[[i]]$loglike))
  }, seq_along(output_raw), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(chain = as.factor(chain))
  
  # initialise final output object
  output_processed <- list(output = df_output)
  
  ## Diagnostics
  output_processed$diagnostics <- list()
  
  # run-times
  run_time <- data.frame(chain = 1:chains,
                         seconds = chain_runtimes)
  output_processed$diagnostics$run_time <- run_time
  
  # Rhat (Gelman-Rubin diagnostic)
  if (chains > 1) {
    rhat_est <- c()
    for (p in 4:19) {
      rhat_est[p - 3] <- df_output %>%
        dplyr::filter(phase == "sampling") %>%
        dplyr::select(c(1, p)) %>%
        gelman_rubin(chains = chains, samples = samples)
    }
    names(rhat_est) <- names(df_output)[4:19]
    output_processed$diagnostics$rhat <- rhat_est
  }
  
  # ESS
  ess_est <- df_output %>%
    dplyr::filter(phase == "sampling") %>%
    dplyr::select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
    apply(2, coda::effectiveSize)
  output_processed$diagnostics$ess <- ess_est
  
  # save output as custom class
  class(output_processed) <- "drjacoby_output"
  
  # return
  return(output_processed)
}

#------------------------------------------------
#' @title Run drjacoby MCMC
#'
#' @description Returns the log-likelihood used within the main MCMC for a fixed
#'   set of parameter values.
#'
#' @param data a list of data elements.
#' @param params a named vector of the free parameters under the model.
#'
#' @export

get_loglike <- function(data,
                        params) {
  
  # get number of weeks in control data and define weekly breaks at hourly
  # resolution
  data_control <- data$data_control
  n_weeks <- ceiling(max(data_control$time.1) / 24 / 7)
  week_breaks <- (0:n_weeks)*7*24
  
  # populate two lists. control_lambda_index tells us which lambda parameter
  # values are needed for each window, and control_lambda_weight tells us the
  # corresponding weighting of these parameters
  control_lambda_index <- list()
  control_lambda_weight <- list()
  for (i in 1:nrow(data_control)) {
    window_days <- data_control$time[i]:(data_control$time.1[i] - 1)
    window_match <- as.numeric(cut(window_days, week_breaks, right = FALSE))
    control_lambda_index[[i]] <- unique(window_match)
    control_lambda_weight[[i]] <- as.vector(table(window_match))
  }
  
  # the eir_adjustment vector is needed in full for the treatment arm, but for
  # the control arm we just need the unique values in this vector and their
  # corresponding weights
  eir_unique <- unique(eir_adjustment)
  eir_weight <- mapply(sum, split(data$ind_weight, f = match(eir_adjustment, eir_unique)))
  
  # data to pass to C++
  args_data <- c(data,
                 list(n_weeks = n_weeks,
                      control_lambda_index = control_lambda_index,
                      control_lambda_weight = control_lambda_weight,
                      eir_unique = eir_unique,
                      eir_weight = eir_weight))
  
  get_loglike_cpp(list(args_data = args_data,
                       params = params))
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

