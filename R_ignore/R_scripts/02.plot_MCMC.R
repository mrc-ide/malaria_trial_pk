# 02.plot_MCMC.R
#
# Author: Bob Verity
# Date: 2023-04-14
#
# Purpose:
# Reads in MCMC output and produces simple plots of model fits against data
#
# ------------------------------------------------------------------

library(tidyverse)
library(pammtools)

# hill function
hill_func <- function(x, min_prob, half_point, hill_power) {
  min_prob + (1 - min_prob) / (1 + (x / half_point)^hill_power)
}

# --------------------------
# read in and process data

# read in MCMC output
mcmc <- readRDS("ignore/outputs/mcmc_raw.rds")

# read in trial data
dat_trial <- read.csv("data/cisse_data.csv")

# read in drug concentration data
dat_drug <- readRDS("data/quadrature_pk.rds")

# get weight of each group*quadrature combination
w_combined <- dat_drug %>%
  dplyr::filter(time == 0) %>%
  mutate(w = weighting * pop_prop) %>%
  pull(w)

# get EIR adjustment (same for all time points)
eir_adjustment <- dat_drug %>%
  dplyr::filter(time == 0) %>%
  pull(eir_adjustment)

# get unique EIRs and weightings
eir_unique <- unique(eir_adjustment)
eir_weight <- mapply(sum, split(w_combined, f = match(eir_adjustment, eir_unique)))

# get drug concentrations into simple matrix
n_ind <- length(w_combined)
drug_mat <- matrix(dat_drug$drug_value.conc, nrow = n_ind, byrow = TRUE)

# drop final column of drug_mat, as this brings the dimension to exactly 13 weeks
drug_mat <- drug_mat[,-ncol(drug_mat)]

n_weeks <- 13
n_hours <- n_weeks * 7 * 24

# --------------------------
# subsample and summarise MCMC output

# subsample MCMC output
n_sub <- 1e2
mcmc_sub <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>%
  sample_n(n_sub)


# for each random sample, calculate prob. uninfected in both control and
# treatment groups
week_vec <- rep(1:n_weeks, each = 7*24)
control_mat <- treat_mat <- matrix(NA, nrow = n_sub, ncol = n_hours)
t0 <- Sys.time()
for (i in 1:n_sub) {
  message(sprintf("%s of %s", i, n_sub))
  
  # control arm
  lambda_vec <- unlist(mcmc_sub[i, sprintf("lambda_%s", 1:n_weeks)])
  r <- cumsum(lambda_vec[week_vec] / 24)
  prob_sus <- 0
  for (j in seq_along(eir_unique)) {
    prob_sus <- prob_sus + eir_weight[j]*exp(-eir_unique[j]*r)
  }
  control_mat[i,] <- 546*(1 - prob_sus)
  
  # treatment arm
  z_mat <- hill_func(drug_mat, mcmc_sub$min_prob[i], mcmc_sub$half_point[i], mcmc_sub$hill_power[i])
  exp_rate <- sweep(z_mat, 2, lambda_vec[week_vec] / 24, "*")
  exp_rate_cumsum <- t(apply(exp_rate, 1, cumsum))
  eir_rate_adjusted <- sweep(exp_rate_cumsum, 1, eir_adjustment, "*")
  treat_mat[i,] <- 542*(1 - colSums(w_combined * exp(-eir_rate_adjusted)))
}
Sys.time() - t0

# summarise into quantiles and get into data.frame
# control mat is a dataframe which for each mcmc sample has the number of infected individuals
# per hour in that arm (given 3 doses). Do this for more samples:

control_q95 <- t(apply(control_mat, 2, quantile_95))
colnames(control_q95) <- sprintf("control_%s", colnames(control_q95))

treat_q95 <- t(apply(treat_mat, 2, quantile_95))
colnames(treat_q95) <- sprintf("treat_%s", colnames(treat_q95))

df_model <- data.frame(time = 0:(nrow(control_q95) - 1)) %>%
  bind_cols(control_q95) %>%
  bind_cols(treat_q95)

saveRDS(df_model, "data/df_model.rds")

# --------------------------
# plot results

# plot against KM data
dat_trial %>%
  group_by(treat_arm) %>%
  dplyr::summarise(time = time.1,
                   n_infected = cumsum(n_infected)) %>%
  mutate(treat_arm = c("control", "treatment")[treat_arm],
         treat_arm = factor(treat_arm)) %>%
  ggplot() + theme_bw() +
  geom_ribbon(aes(x = time / 24, ymin = control_Q2.5, ymax = control_Q97.5),
              alpha = 0.2, fill = "red", data = df_model) +
  geom_line(aes(x = time / 24, y = control_Q50), alpha = 0.2, col = "red", data = df_model) +
  geom_ribbon(aes(x = time / 24, ymin = treat_Q2.5, ymax = treat_Q97.5),
              alpha = 0.2, fill = "blue", data = df_model) +
  geom_line(aes(x = time / 24, y = treat_Q50), alpha = 0.2, col = "blue", data = df_model) +
  #geom_step(aes(x = time / 24, y = n_infected, col = treat_arm), direction = "hv") +
  geom_point(aes(x = time / 24, y = n_infected, col = treat_arm), size = 0.8) +
  xlab("Time (days)") + ylab("Number infected")

## plot data vs fitting
df_model_long <- df_model %>% tidyr::pivot_longer(control_Q2.5:treat_Q97.5,
                                                  names_to = "arm_Q", 
                                                  values_to = "cum_infected") %>% 
  dplyr::mutate(arm = stringr::str_split_fixed(arm_Q, pattern = "_", n = 2)[,1],
                quantile = stringr::str_split_fixed(arm_Q, pattern = "_", n = 2)[,2]) %>%
  dplyr::select(time, arm, quantile, cum_infected) %>%
  dplyr::arrange(arm, quantile, time) %>%
  dplyr::mutate(arm = if_else(arm == "control", 1, 2))

df_model_50 <- df_model_long %>% 
  dplyr::filter(quantile == "Q50")

compare_mcmc_arm <- function(real_data = dat_trial, drug_arm, q = "Q50", mcmc_output = df_model_long) {
  #' input: real data to gain time steps, which treatment arm, output from the mcmc (long format with all arms and qs)
  #' process: extract time steps, initialise the number of patients, convert mcmc data to infections in time intervals from data
  #' output: dataframe comparible to the original trial data (n_infected/n_patients)
  treat_data <- real_data %>%
    dplyr::filter(treat_arm == drug_arm)
  max_time <- treat_data$time.1[nrow(treat_data)]
  times <- c(treat_data$time, max_time)
  N_patients <- treat_data$n_patients[1]
  N_remaining <- N_patients
  
  mcmc <- mcmc_output %>%
    dplyr::filter(arm == drug_arm) %>%
    dplyr::filter(quantile == q)
  
  df_mcmc <- treat_data
  new_infections <- 0
  df_mcmc$n_patients <- NA
  df_mcmc$n_infected <- NA
  
  for (i in 1:nrow(df_mcmc)) {
    N_remaining <- N_remaining - new_infections
    time0 <- df_mcmc$time[i]
    time1 <- df_mcmc$time.1[i]
    new_infections <- mcmc$cum_infected[mcmc$time == time1] - mcmc$cum_infected[mcmc$time == time0]
    ## simulated df
    df_mcmc$n_infected[i] <- mcmc$cum_infected[mcmc$time == time1] - mcmc$cum_infected[mcmc$time == time0]
    df_mcmc$n_patients[i] <- N_remaining
    df_mcmc$treat_arm <- as.factor(drug_arm)
  }
  
  df_mcmc <- df_mcmc %>%
    dplyr::select(c(n_patients, n_infected, time, time.1, treat_arm))
  return(df_mcmc)
}

df_mcmc_control <- compare_mcmc_arm(drug_arm = 1)
df_mcmc_treat <- compare_mcmc_arm(drug_arm = 2)
df_mcmc <- rbind(df_mcmc_control, df_mcmc_treat) %>%
  dplyr::mutate(method = "mcmc") %>%
  dplyr::mutate(compare = n_infected/n_patients)

dat_trial <- dat_trial %>% 
  dplyr::select(names(df_mcmc_control)) %>%
  dplyr::mutate(method = "data") %>%
  dplyr::mutate(compare = n_infected/n_patients)

rbind(df_mcmc, dat_trial) %>%
  dplyr::select(time.1:compare) %>%
  tidyr::pivot_wider(names_from = method, values_from = compare) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = data, y = mcmc, col = treat_arm)) +
  geom_abline(slope = 1, intercept = 0, lty = 2)

##################################################################################################
## NEW PLOTTING METHODOLOGY ##
head(df_model_long)

df_model_long <- df_model_long %>%
  dplyr::group_by(quantile, arm) %>% 
  dplyr::mutate(incidence = cum_infected - lag(cum_infected))

df_long_1 <- df_model_long %>%
  dplyr::filter(arm == 1)
df_long_2 <- df_model_long %>%
  dplyr::filter(arm == 2)

## from the trial data, we want to sum up the incidence within the intervals of the KM data
dat_trial <- read.csv("data/cisse_data.csv")
head(dat_trial)

dat_trial_1 <- dat_trial %>%
  dplyr::filter(treat_arm == 1)
dat_trial_2 <- dat_trial %>%
  dplyr::filter(treat_arm == 2)

time_1 <- c(dat_trial_1$time,dat_trial_1$time.1[nrow(dat_trial_1)])
time_2 <- c(dat_trial_2$time,dat_trial_2$time.1[nrow(dat_trial_2)])

## arm 1
df_sim_1 <- dat_trial_1 %>%
  dplyr::select(time, time.1, treat_arm)

for(i in 1:nrow(df_sim_1)) {
  df <- df_long_1 %>%
    dplyr::filter(quantile == "Q50") %>%
    dplyr::filter(time < df_sim_1$time.1[i]) %>%
    dplyr::filter(time >= df_sim_1$time[i]) 
  
  df_sim_1$incidence_50[i] <- sum(df$incidence, na.rm = TRUE)
}

for(i in 1:nrow(df_sim_1)) {
  df <- df_long_1 %>%
    dplyr::filter(quantile == "Q2.5") %>%
    dplyr::filter(time < df_sim_1$time.1[i]) %>%
    dplyr::filter(time >= df_sim_1$time[i]) 
  
  df_sim_1$incidence_2.5[i] <- sum(df$incidence, na.rm = TRUE)
}

for(i in 1:nrow(df_sim_1)) {
  df <- df_long_1 %>%
    dplyr::filter(quantile == "Q97.5") %>%
    dplyr::filter(time < df_sim_1$time.1[i]) %>%
    dplyr::filter(time >= df_sim_1$time[i]) 
  
  df_sim_1$incidence_97.5[i] <- sum(df$incidence, na.rm = TRUE)
}
df_sim_1$data <- dat_trial_1$n_infected
df_sim_1 <- df_sim_1 %>% dplyr::select(time, time.1, treat_arm, incidence_2.5, incidence_50, incidence_97.5, data)

ggplot(data = df_sim_1) + theme_bw() + 
  geom_step(aes(x = time, y = data)) +
  geom_step(aes(x = time, y = incidence_50), linetype = "dashed", colour = "blue") + 
  geom_stepribbon(aes(x = time, ymin = incidence_2.5, ymax = incidence_97.5), fill = "blue", alpha = 0.2)

## arm 2
df_sim_2 <- dat_trial_2 %>%
  dplyr::select(time, time.1, treat_arm)

for(i in 1:nrow(df_sim_2)) {
  df <- df_long_2 %>%
    dplyr::filter(quantile == "Q50") %>%
    dplyr::filter(time < df_sim_2$time.1[i]) %>%
    dplyr::filter(time >= df_sim_2$time[i]) 
  
  df_sim_2$incidence_50[i] <- sum(df$incidence, na.rm = TRUE)
}

for(i in 1:nrow(df_sim_2)) {
  df <- df_long_2 %>%
    dplyr::filter(quantile == "Q2.5") %>%
    dplyr::filter(time < df_sim_2$time.1[i]) %>%
    dplyr::filter(time >= df_sim_2$time[i]) 
  
  df_sim_2$incidence_2.5[i] <- sum(df$incidence, na.rm = TRUE)
}

for(i in 1:nrow(df_sim_2)) {
  df <- df_long_2 %>%
    dplyr::filter(quantile == "Q97.5") %>%
    dplyr::filter(time < df_sim_2$time.1[i]) %>%
    dplyr::filter(time >= df_sim_2$time[i]) 
  
  df_sim_2$incidence_97.5[i] <- sum(df$incidence, na.rm = TRUE)
}
df_sim_2$data <- dat_trial_2$n_infected
df_sim_2 <- df_sim_2 %>% dplyr::select(time, time.1, treat_arm, incidence_2.5, incidence_50, incidence_97.5, data)

ggplot(data = df_sim_2) + theme_bw() + 
  geom_step(aes(x = time, y = data)) +
  geom_step(aes(x = time, y = incidence_50), linetype = "dashed", colour = "blue") + 
  geom_stepribbon(aes(x = time, ymin = incidence_2.5, ymax = incidence_97.5), fill = "blue", alpha = 0.2)

df_sim <- rbind(df_sim_1, df_sim_2)

ggplot(data = df_sim) + theme_bw() + 
  geom_step(aes(x = time, y = data)) +
  geom_step(aes(x = time, y = incidence_50), linetype = "dashed", colour = "blue") + 
  geom_stepribbon(aes(x = time, ymin = incidence_2.5, ymax = incidence_97.5), fill = "blue", alpha = 0.2) + 
  facet_grid(treat_arm ~ .)
