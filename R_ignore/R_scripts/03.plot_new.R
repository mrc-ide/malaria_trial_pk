## plotting script from scratch to try and fix the issues with the current plotting script
library(tidyverse)
library(pammtools)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra) 

# hill function
hill_func <- function(x, min_prob, half_point, hill_power) {
  min_prob + (1 - min_prob) / (1 + (x / half_point)^hill_power)
}

# --------------------------
# read in and process data

# read in MCMC output
mcmc <- readRDS("ignore/outputs/mcmc_raw_final.rds")

# read in trial data
dat_trial <- read.csv("data/cisse_data.csv") #%>%
  # dplyr::filter(time.1 < 2000) ## appended data

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
drug_mat <- drug_mat[,-ncol(drug_mat)] # eat column is an hour, each row is a quadrature group

n_weeks <- 13
n_hours <- n_weeks * 7 * 24
## appended data only
# n_hours <- max(dat_trial$time.1)
# n_weeks <- n_hours/(24*7)

# --------------------------
# subsample and summarise MCMC output

# subsample MCMC output
n_sub <- 1e3
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
  prob_sus <- prob_sus[1:ncol(control_mat)]
  control_mat[i,] <- 546*(1 - prob_sus)
  
  # treatment arm
  z_mat <- hill_func(drug_mat, mcmc_sub$min_prob[i], mcmc_sub$half_point[i], mcmc_sub$hill_power[i])
  exp_rate <- sweep(z_mat, 2, lambda_vec[week_vec] / 24, "*")
  exp_rate_cumsum <- t(apply(exp_rate, 1, cumsum))
  eir_rate_adjusted <- sweep(exp_rate_cumsum, 1, eir_adjustment, "*")[,1:ncol(treat_mat)]
  treat_mat[i,] <- 542*(1 - colSums(w_combined * exp(-eir_rate_adjusted)))
}
Sys.time() - t0

## look at data -- cumulative infections in each of the sub-sample (increase n_sub for final run)
# View(control_mat)
# View(treat_mat)

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

## bug: runs line by line but not in function. using line by line version for first pass:
control_summary <- quantiles_mcmc_arm(drug_arm = 1, mcmc_sample = control_mat)
treat_summary <- quantiles_mcmc_arm(drug_arm = 2, mcmc_sample = treat_mat)

saveRDS(control_summary, "control_summary.rds")
saveRDS(treat_summary, "treat_summary.rds")

mcmc_summary <- rbind(control_summary, treat_summary)
saveRDS(mcmc_summary, "mcmc_summary.rds")

# convert data to rates
df_trial <- dat_trial %>%
  dplyr::group_by(treat_arm) %>%
  dplyr::mutate(rate = n_infected/n_patients/(time.1-time))
df_trial$treat_arm <- as.factor(df_trial$treat_arm)

arm_labs <- c("Placebo", "SP+AS")
names(arm_labs) <- c(1,2)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col_scale_vals <- gg_color_hue(n = 2)
names(col_scale_vals) <- as.factor(c(1,2))
col_scale <- scale_colour_manual(name = "treat_arm", values = col_scale_vals)
fill_scale <- scale_fill_manual(name = "treat_arm", values = col_scale_vals)

mcmc_summary_treat <- mcmc_summary %>% 
  dplyr::filter(treat_arm == 2) %>%
  dplyr::mutate(facet_factor = "Treatment")
mcmc_summary_combined <- mcmc_summary %>% 
  dplyr::mutate(facet_factor = "Combined")
mcmc_summary_final <- rbind(mcmc_summary_treat, mcmc_summary_combined)

df_treat <- df_trial %>% 
  dplyr::filter(treat_arm == 2) %>%
  dplyr::mutate(facet_factor = "Treatment")
df_combined <- df_trial %>% 
  dplyr::mutate(facet_factor = "Combined")
df_final <- rbind(df_treat, df_combined)

ggplot() + geom_step(data = df_final, aes(x = time/24, y = rate * 1000, col = treat_arm)) + 
  facet_grid(facet_factor ~ ., scales = "free_y") +
  theme_bw() + theme(strip.text.y = element_blank()) +
  geom_stepribbon(data = mcmc_summary_final, 
                  aes(x = time/24, ymin = rate_2.5* 1000, ymax = rate_97.5* 1000, fill = treat_arm), alpha = 0.2) + 
  geom_step(data = mcmc_summary_final, aes(x = time/24, y = rate_50* 1000, col = treat_arm), linetype = 2) + 
  scale_color_discrete(labels = arm_labs) + 
  scale_fill_discrete(labels = arm_labs) +
  guides(fill = guide_legend(title = "Treatment arm"),
         col = guide_legend(title = "Treatment arm")) +
  labs(x = "Time (days)", y = "Infection rate per 1,000 children at risk per day ") + 
  geom_vline(xintercept = c(0, 28, 56), linetype = "dotted")

ggsave("output/mcmc_fit.png", dpi = 500, width = 20, height = 14, units = "cm")
