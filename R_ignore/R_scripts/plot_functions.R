# hill function
hill_func <- function(x, min_prob, half_point, hill_power) {
  ret <- min_prob + (1 - min_prob) / (1 + (x / half_point)^hill_power)
  ret
}

## convert matrices into correct format before finding quantiles
quantiles_mcmc_arm <- function(real_data = dat_trial, drug_arm, mcmc_sample = arm_mat) {
  #' input: real data to gain time steps, which treatment arm, output from the mcmc (wide format with cum inf each hour per sampled params)
  #' process: extract time steps, convert mcmc data to infections in time intervals from data
  #' output: dataframe comparible to the original trial data (n_infected/n_patients) with quantiles
  treat_data <- real_data %>%
    dplyr::filter(treat_arm == drug_arm)
  max_time <- treat_data$time.1[nrow(treat_data)]
  times <- c(treat_data$time, max_time)
  N_patients <- treat_data$n_patients[1]
  N_remaining <- N_patients
  
  df_mcmc <- treat_data[rep(seq_len(nrow(treat_data)), n_sub), ]
  df_mcmc$sample <- rep(1:n_sub, each = nrow(treat_data))
  df_mcmc$n_patients <- as.numeric(0)
  df_mcmc$n_infected <- 0
  
  
  df <- df_mcmc %>% 
    dplyr::group_by(sample) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n_infected = mcmc_sample[sample,time.1] - mcmc_sample[sample,time+1],
                  n_patients = N_patients - mcmc_sample[sample,time.1]) %>%
    dplyr::mutate(rate = n_infected/n_patients/(time.1-time))
  
  df_summary <- df %>% 
    dplyr::group_by(time) %>%
    dplyr::summarise(n_inf_2.5 = quantile(n_infected, 0.025, na.rm = TRUE),
                     n_inf_50 = quantile(n_infected, 0.5, na.rm = TRUE),
                     n_inf_97.5 = quantile(n_infected, 0.975, na.rm = TRUE),
                     n_pat_2.5 = quantile(n_patients, 0.025, na.rm = TRUE),
                     n_pat_50 = quantile(n_patients, 0.5, na.rm = TRUE),
                     n_pat_97.5 = quantile(n_patients, 0.975, na.rm = TRUE),
                     rate_2.5 = quantile(rate, 0.025, na.rm = TRUE),
                     rate_50 = quantile(rate, 0.5, na.rm = TRUE),
                     rate_97.5 = quantile(rate, 0.975, na.rm = TRUE))
  df_summary$time.1 <- treat_data$time.1
  df_summary$treat_arm <- as.factor(drug_arm)
  
  df_summary <- df_summary %>% 
    dplyr::relocate(time.1, .before = n_inf_2.5)
  
  return(df_summary)
  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

age_to_group <- function(age) {
  if (age < 12) {
    group <- "0 - 1 years"
  }

  if (age >= 12 & age < 24) {
    group <- "1 - 2 years"
  }

  if (age >= 24 & age < 36) {
    group <- "2 - 3 years"
  }

  if (age >= 36 & age < 48) {
    group <- "3 - 4 years"
  }

  if (age >= 48 & age < 60) {
    group <- "4 - 5 years"
  }
  group
}


group_to_zscore <- function(group) {
  if (group %% 2 == 0) {
    zscore = "not malnourished"
  }

  if (group %% 2 == 1) {
    zscore = "malnourished"
  }

  zscore
}

# pharmacodynamics
median_pd <- function(conc_vec, mcmc_cri) {
  eff_vec <- 1 - hill_func(x = conc_vec,
                           min_prob = mcmc_cri$median[1],
                           half_point = mcmc_cri$median[2],
                           hill_power = mcmc_cri$median[3])
  
  cri_df <- data.frame(concentration = conc_vec, 
                       efficacy = eff_vec)
  
  return(cri_df)
}

sample_pd <- function(mcmc_output, num_samples = 100, conc_vec) {
  mcmc_post <- mcmc_output %>% 
    dplyr::filter(phase == "sampling") %>%
    dplyr::slice_sample(n = num_samples) %>% # select 100 random iterations from sampling phase 
    dplyr::select(min_prob:hill_power)
  
  df <- data.frame(time = numeric(0), concentration = numeric(0),
                   efficacy = numeric(0), sample = numeric(0)) 
  len <- length(conc_vec)
  
  for(i in 1:nrow(mcmc_post)){ 
    y <- data.frame(concentration = numeric(len),
                    efficacy = numeric(len), sample = numeric(len))
    y$concentration <- conc_vec
    y$efficacy <- 1 - hill_func(x = y$concentration,
                                min_prob = mcmc_post$min_prob[i],
                                half_point = mcmc_post$half_point[i],
                                hill_power = mcmc_post$hill_power[i])
    
    y$sample <- i
    
    df <- rbind(df, y)
  }
  return(df)
}

# pharmacokinetics
median_pk <- function(quad_pk) {
  median_df <- quad_pk %>%
    dplyr::group_by(time) %>%
    dplyr::reframe(median_conc = reldist::wtd.quantile(drug_value.conc, q = 0.5,
                                                       weight = weighting*pop_prop)) %>%
    dplyr::distinct()
  
  max_data <- max(median_df$median_conc)
  ind <- which(median_df$median_conc == max_data)
  data_subset <- median_df[ind:nrow(median_df),]
  data_subset[,2] <- log(data_subset[,2])
  names(data_subset)[2] <- "log_concentration"
  
  fit <- lm(formula = log_concentration ~ time, data = data_subset)
  
  c_0 <- exp(fit$coefficients[1]) 
  k <- -(fit$coefficients[2]) 
  
  df <- median_df
  df$median_conc <- c_0 * exp(-k * df$time)
  
  return(df)
}

sample_posterior <- function(mcmc_output, num_samples = 100, median_df) {
  mcmc_post <- mcmc_output %>% 
    dplyr::filter(phase == "sampling") %>%
    dplyr::slice_sample(n = num_samples) %>% # select 100 random iterations from sampling phase 
    dplyr::select(min_prob:hill_power)

  df <- data.frame(time = numeric(0), concentration = numeric(0),
                   efficacy = numeric(0), sample = numeric(0)) 
  len <- nrow(median_df)
  
  for(i in 1:nrow(mcmc_post)){ 
    y <- data.frame(time = numeric(len), concentration = numeric(len),
                    efficacy = numeric(len), sample = numeric(len))
    y$time <- median_df$time
    y$concentration <- median_df$median_conc
      
    for(j in 1:nrow(y)) {
      y$efficacy[j] <- 1 - hill_func(x = y$concentration[j]/1000,
                                     min_prob = mcmc_post$min_prob[i],
                                     half_point = mcmc_post$half_point[i],
                                     hill_power = mcmc_post$hill_power[i])
    }
    y$sample <- i
    
    df <- rbind(df, y)
  }
  return(df)
}

median_posterior <- function(mcmc_summary, median_df) { 
  min_prob <- mcmc_summary$median[1]
  half_point <- mcmc_summary$median[2]
  hill_power <- mcmc_summary$median[3]
  
  df <- median_df
  df$efficacy <- numeric(nrow(df))
  for(i in 1:nrow(df)) {
    df$efficacy[i] <- 1 - hill_func(x = df$median_conc[i]/1000,
                                min_prob = min_prob,
                                half_point = half_point,
                                hill_power = hill_power)
  }
  return(df)
}

sample_pe <- function(mcmc_output, num_samples, median_pk) {
  
  mcmc_post <- mcmc_output %>% 
    dplyr::filter(phase == "sampling") %>%
    dplyr::slice_sample(n = num_samples) %>% # select 100 random iterations from sampling phase 
    dplyr::select(min_prob:hill_power)
  
  df <- data.frame(time = numeric(0),
                   concentration = numeric(0),
                   sample = numeric(0))
  for(i in 1:num_samples) {
    y <- data.frame(time = numeric(nrow(median_pk)),
                    concentration = numeric(nrow(median_pk)),
                    sample = numeric(nrow(median_pk)))
    y$time <- median_pk$time
    y$concentration <- median_pk$median_conc
    y$efficacy <- 1 - hill_func(x = y$concentration,
                                min_prob = mcmc_post$min_prob[i],
                                half_point = mcmc_post$half_point[i],
                                hill_power = mcmc_post$hill_power[i])
    y$sample <- i
    df <- rbind(df, y)
  }
  return(df)
  
}
