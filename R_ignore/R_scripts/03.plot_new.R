## plotting script from scratch to try and fix the issues with the current plotting script
library(tidyverse)
library(pammtools)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra) 
library(reldist)

source("R_ignore/R_scripts/plot_functions.R")

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

quadrature_pk_1 <- readRDS("data/quadrature_pk_1.rds")

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
n_sub <- 1e1
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

# control_summary <- quantiles_mcmc_arm(drug_arm = 1, mcmc_sample = control_mat)
# treat_summary <- quantiles_mcmc_arm(drug_arm = 2, mcmc_sample = treat_mat)
# 
# saveRDS(control_summary, "control_summary.rds")
# saveRDS(treat_summary, "treat_summary.rds")
# 
# mcmc_summary <- rbind(control_summary, treat_summary)
# saveRDS(mcmc_summary, "mcmc_summary.rds")
mcmc_summary <- readRDS("mcmc_summary.rds")

# convert data to rates
df_trial <- dat_trial %>%
  dplyr::group_by(treat_arm) %>%
  dplyr::mutate(rate = n_infected/n_patients/(time.1-time))
df_trial$treat_arm <- as.factor(df_trial$treat_arm)

arm_labs <- c("Placebo", "SP+AS")
names(arm_labs) <- c(1,2)

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

mcmc_cri <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>% 
  tidyr::pivot_longer(cols = lambda_1:hill_power, names_to = "parameter", values_to = "value") %>% 
  dplyr::group_by(parameter)  %>%
  dplyr::mutate(scaling = if_else(grepl("lambda", parameter), 1000, 1)) %>%
  dplyr::reframe(scaling = scaling, 
                 median = median(value) * scaling,
                 lower_cri = quantile(value, probs = 0.025) * scaling,
                 upper_cri = quantile(value, probs = 0.975) * scaling) %>%
  distinct() %>%
  dplyr::arrange(factor(parameter, 
                        levels = c("min_prob", "half_point", "hill_power",
                                   paste0("lambda_", 1:13))))

## diagnostics
mcmc$diagnostics

## figure 3 with the 6 panels, based upon the pharmacokinetics by age
quadrature_pk <- readRDS("data/quadrature_pk_1.rds")
concentration_vec <- seq(0, 50, by = 0.1)
## pharamacodynamics
pd_median <- median_pd(conc_vec = concentration_vec,
                       mcmc_cri = mcmc_cri)
pd_cri <- sample_pd(mcmc_output = mcmc$output, num_samples = 1000,
                    conc_vec = concentration_vec)
pd_plot <- ggplot() +
  geom_line(data = pd_cri, aes(x = concentration, y = efficacy, group = sample),
              alpha = 0.1, col = "darkgrey") +
  geom_line(data = pd_median, aes(x = concentration, y = efficacy), lwd = 1) +
  scale_x_log10() + theme_bw() +
  labs(x = "log concentration (\u03bcl / ml)", y = "protective efficacy")

## protective efficacy
median_df <- median_pk(quad_pk)
mcmc_post <- sample_posterior(mcmc_output = mcmc$output, num_samples = 1000, 
                              median_df = median_df)
head(mcmc_post)

mcmc_median <- median_posterior(mcmc_summary = mcmc_cri, median_df = median_df)

pk_plot <- ggplot() + 
  geom_line(data = mcmc_post, aes(x = time/24, y = efficacy, group = sample),
            alpha = 0.05, col = "darkgrey") +
  theme_bw() +
  geom_line(data = mcmc_median, aes(x = time/24, y = efficacy), lwd = 1) +
  xlim(c(0.5/24, 60)) + labs(x = "time (days)") 


quadrature_pk_age_nut <- quad_pk %>%
  dplyr::group_by(individual) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(age_group = age_to_group(age),
                z_score = group_to_zscore(as.numeric(group)))

pk_age_med <- quadrature_pk_age_nut %>% 
  dplyr::group_by(time, age_group) %>% 
  dplyr::summarise(lower_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.025),
                   median = median(drug_value.conc),
                   upper_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.975))
pk_nut_med <- quadrature_pk_age_nut %>% 
  dplyr::group_by(time, z_score) %>% 
  dplyr::summarise(lower_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.025),
                   median = median(drug_value.conc),
                   upper_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.975))

quadrature_pk_cri <- quad_pk %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(lower_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.025),
                   median = median(drug_value.conc),
                   upper_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.975))

ggplot(data = quadrature_pk_cri, aes(x = time/24, y = median)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = lower_cri, ymax = upper_cri), fill = "blue", alpha = .1) +
  theme_bw() + xlim(c(0, 50)) + 
  labs(x = "days", y = "median sulfadoxine concentration (\u03bcg / ml)")

pk_age_med$age_group <- factor(pk_age_med$age_group)
pk_nut_med$z_score <- factor(pk_nut_med$z_score)

ggplot() + 
  geom_line(data = pk_age_med, aes(x = time/24, y = median, col = age_group)) + 
  theme_bw() + xlim(c(0,50)) + 
  geom_line(data = quadrature_pk_cri, aes(x = time/24, y = median), 
            col = "black", linetype = 2) + 
  geom_ribbon(data = quadrature_pk_cri, aes(x = time/24, y = median, 
                                            ymin = lower_cri, ymax = upper_cri),
              alpha = .1) +
  labs(x = "days", y = "median sulfadoxine concentration (\u03bcg / ml)")

ggplot() + 
  geom_line(data = pk_nut_med, aes(x = time/24, y = median, col = z_score)) + 
  theme_bw() + xlim(c(0,50)) + 
  geom_line(data = quadrature_pk_cri, aes(x = time/24, y = median), 
            col = "black", linetype = 2) + 
  geom_ribbon(data = quadrature_pk_cri, aes(x = time/24, y = median, 
                                            ymin = lower_cri, ymax = upper_cri),
              alpha = .1) +
  labs(x = "days", y = "median sulfadoxine concentration (\u03bcg / ml)")
