# 02.plot_MCMC.R
#
# Author: Bob Verity
# Date: 2023-04-14
#
# Purpose:
# Reads in MCMC output and produces simple plots of model fits against data
#
# ------------------------------------------------------------------

# hill function
hill_func <- function(x, min_prob, half_point, hill_power) {
  min_prob + (1 - min_prob) / (1 + (x / half_point)^hill_power)
}

# --------------------------
# read in and process data

# read in MCMC output
mcmc <- readRDS("R_ignore/outputs/mcmc_raw.rds")

# read in trial data
dat_trial <- read.csv("R_ignore/data/cisse_data.csv")

# read in drug concentration data
dat_drug <- readRDS("R_ignore/data/quadrature_pk.rds")

# get weight of each group*quadrature combination
w_combined <- dat_drug %>%
  dplyr::filter(time == 0) %>%
  mutate(w = weighting * pop_prop) %>%
  pull(w)
w_combined <- w_combined / sum(w_combined)

# get drug concentrations into simple matrix
n_ind <- length(w_combined)
drug_mat <- matrix(dat_drug$drug_value.conc, nrow = n_ind, byrow = TRUE)

# --------------------------
# subsample and summarise MCMC output

# subsample MCMC output
n_sub <- 1e2
mcmc_sub <- mcmc %>%
  dplyr::filter(phase == "sampling") %>%
  sample_n(n_sub)

# for each random sample, calculate prob. uninfected in both control and
# treatment groups
t_vec <- min(dat_drug$time):max(dat_drug$time)
control_mat <- treat_mat <- matrix(NA, nrow = n_sub, ncol = length(t_vec))
for (i in 1:n_sub) {
  message(sprintf("%s of %s", i, n_sub))
  
  # control arm
  control_mat[i,] <- 546*(1 - exp(-mcmc_sub$lambda[i] / 24 * t_vec))
  
  # treatment arm
  z_mat <- hill_func(drug_mat, mcmc_sub$min_prob[i], mcmc_sub$half_point[i], mcmc_sub$hill_power[i])
  z_cumsum <- t(apply(z_mat, 1, cumsum))
  treat_mat[i,] <- 542*(1 - colSums(w_combined * exp(-mcmc_sub$lambda[i] / 24 * z_cumsum)))
}

# summarise into quantiles and get into data.frame
control_q95 <- t(apply(control_mat, 2, quantile_95))
colnames(control_q95) <- sprintf("control_%s", colnames(control_q95))

treat_q95 <- t(apply(treat_mat, 2, quantile_95))
colnames(treat_q95) <- sprintf("treat_%s", colnames(treat_q95))

df_model <- data.frame(time = t_vec) %>%
  bind_cols(control_q95) %>%
  bind_cols(treat_q95)

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
  geom_step(aes(x = time / 24, y = n_infected, col = treat_arm), direction = "hv") +
  geom_point(aes(x = time / 24, y = n_infected, col = treat_arm), size = 0.8) +
  xlab("Time (days)") + ylab("Number infected")



