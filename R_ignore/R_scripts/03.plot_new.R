## plotting script from scratch to try and fix the issues with the current plotting script
library(tidyverse)
library(pammtools)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra) 
library(reldist)

source("R_ignore/R_scripts/plot_functions.R")
devtools::load_all(".")
old <- ggplot2::theme_set(theme_bw(base_size = 12)) 

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

control_mat[1:10,1:10]


control_summary <- quantiles_mcmc_arm(drug_arm = 1, mcmc_sample = control_mat)
treat_summary <- quantiles_mcmc_arm(drug_arm = 2, mcmc_sample = treat_mat)
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
  labs(x = "log concentration (\u03bcg / ml)", y = "Protective efficacy")

## protective efficacy
median_df <- median_pk(quadrature_pk_1)
mcmc_post <- sample_posterior(mcmc_output = filter(mcmc$output, phase == "sampling"), num_samples = 1000, 
                              median_df = median_df)
head(mcmc_post)

mcmc_median <- median_posterior(mcmc_summary = mcmc_cri, median_df = median_df)

pe_plot <- ggplot() + 
  geom_line(data = mcmc_post, aes(x = time/24, y = efficacy, group = sample),
            alpha = 0.05, col = "darkgrey") +
  theme_bw() +
  geom_line(data = mcmc_median, aes(x = time/24, y = efficacy), lwd = 1) +
  xlim(c(0.5/24, 60)) + labs(x = "Time (days)", y = "Concentration (\u03bcg / ml)") 


quadrature_pk_age_nut <- quadrature_pk_1 %>%
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

quadrature_pk_cri <- quadrature_pk_1 %>%
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

pk_age_med <- pk_age_med %>% 
  dplyr::group_by(age_group) %>%
  dplyr::mutate(efficacy = 1 - hill_func(x = median/1000,
                                         min_prob = mcmc_cri$median[1],
                                         half_point = mcmc_cri$median[2],
                                         hill_power = mcmc_cri$median[3]))
pk_nut_med <- pk_nut_med %>% 
  dplyr::group_by(z_score) %>%
  dplyr::mutate(efficacy = 1 - hill_func(x = median/1000,
                                         min_prob = mcmc_cri$median[1],
                                         half_point = mcmc_cri$median[2],
                                         hill_power = mcmc_cri$median[3]))

quad_df <- quadrature_pk_1 %>%
  dplyr::select(c(individual, time, drug_value.conc, weighting, pop_prop))

mcmc_output <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>%
  dplyr::select(c(min_prob, half_point, hill_power)) 

# profile_pe <- profvis::profvis(pe_cri(mcmc_output[1:10,], quadrature_pk_1))
pe_sample <- pe_cri(mcmc_output, quad_df)
pe_cri <- pe_sample %>% 
  dplyr::group_by(time) %>%
  dplyr::reframe(lower_cri = quantile(efficacy, 0.025),
                 median = median(efficacy),
                 upper_cri = quantile(efficacy, 0.975))
ggplot(pe_cri) + geom_line(aes(x=time/24, y = median))+
  geom_ribbon(aes(x = time/24, ymin = lower_cri, ymax = upper_cri))

saveRDS(pe_sample, "data/pe_sample.rds")
saveRDS(pe_cri, "data/pe_cri.rds")

ggplot() + geom_line(data = pe_sample, aes(x = time/24, y = efficacy, group = sample), 
                     alpha = 0.15, col = "darkgrey") +
  ylim(c(0,1)) + theme_bw() + xlim(c(0.5, 60)) +
  geom_line(data = pk_age_med, aes(x = time/24, y = efficacy, col = age_group), linewidth = 1) + 
  labs(x = "time (days)", y = "protective efficacy") +
  guides(col = guide_legend(title = "Age group")) 

ggplot() + geom_line(data = pe_sample, aes(x = time/24, y = efficacy, group = sample), 
                     alpha = 0.15, col = "darkgrey") +
  ylim(c(0,1)) + theme_bw() + xlim(c(0.5, 60)) +
  geom_line(data = pk_nut_med, aes(x = time/24, y = efficacy, col = z_score), linewidth = 1) + 
  labs(x = "time (days)", y = "protective efficacy") +
  guides(col = guide_legend(title = "Nutrition status"))

# drug separation 
aq_wide <- aq_df %>%
  tidyr::pivot_wider(names_from = val, values_from = efficacy)

# time in days
ggplot(aq_wide, aes(x = time)) + geom_line(aes(y = med)) +
  geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.2, fill = "blue")

# time in hours
sp_pe_cri <- pe_sample %>%
  dplyr::group_by(time) %>%
  dplyr::reframe(lower_cri = quantile(efficacy, probs = 0.025),
                 median = median(efficacy),
                 upper_cri = quantile(efficacy, probs = 0.975)) %>%
  dplyr::mutate(time = time/24)

ggplot(data = sp_pe_cri) + geom_ribbon(aes(x = time, 
                                           ymin = lower_cri, ymax = upper_cri), 
                                       fill = "red", alpha = 0.2) +
  geom_line(aes(x = time, y = median), lwd = 0.4) +
  labs(x = "Time (days)", y = "Protective efficacy") + xlim(c(0.5, 60))

names(aq_wide) <- c("time", "low_aq", "med_aq", "upp_aq")
names(sp_pe_cri) <- c("time", "low_sp", "med_sp", "upp_sp")

drug_sep <- dplyr::inner_join(aq_wide, sp_pe_cri) %>%
  dplyr::mutate(med_spaq = 1 - ((1-med_aq)*(1-med_sp)),
                low_spaq = 1 - ((1-low_aq)*(1-low_sp)),
                upp_spaq = 1 - ((1-upp_aq)*(1-upp_sp)))
drug_sep$resistance <- "Triple"

## don't think I can just take the lower and upper CI's on the parameters to compute the CI on the distribution

sp_res <- read.csv("data/quint_cri.csv")

names(aq_wide) <- c("time", "low_aq", "med_aq", "upp_aq")
names(sp_res) <- c("time", "low_sp", "med_sp", "upp_sp")
drug_sep_res <- dplyr::inner_join(aq_wide, sp_res) %>%
  dplyr::mutate(med_spaq = 1 - ((1-med_aq)*(1-med_sp)),
                low_spaq = 1 - ((1-low_aq)*(1-low_sp)),
                upp_spaq = 1 - ((1-upp_aq)*(1-upp_sp)))

ggplot(drug_sep_res) + 
  geom_line(aes(x = time, y = med_aq), col = "blue") +
  geom_ribbon(aes(x = time, ymin = low_aq, ymax = upp_aq), alpha = 0.2, fill = "blue") +
  geom_line(aes(x = time, y = med_sp), col = "red")+
  geom_ribbon(aes(x = time, ymin = low_sp, ymax = upp_sp), alpha = 0.2, fill = "red") +
  geom_line(aes(x = time, y = med_spaq)) +
  geom_ribbon(aes(x = time, ymin = low_spaq, ymax = upp_spaq), alpha = 0.5, fill = "darkgrey") + 
  geom_point(x = 28, y = 0.90) +
  geom_point(x = 28, y = 0.83) +
  geom_point(x = 28, y = 0.77) +
  labs(x = "Time (days)", y = "Protective efficacy") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))

drug_sep_res$resistance <- "Quintuple"

drug_sep_long <- drug_sep %>%
  tidyr::pivot_longer(low_aq:upp_spaq, names_to = "int_drug", values_to = "efficacy") %>%
  dplyr::filter(time > 0)
drug_res_long <- drug_sep_res %>%
  tidyr::pivot_longer(low_aq:upp_spaq, names_to = "int_drug", values_to = "efficacy") %>%
  dplyr::filter(time > 0)
drug_res <- rbind(drug_sep_long, drug_res_long) %>%
  rowwise() %>%
  dplyr::mutate(interval = unlist(strsplit(int_drug, "_"))[1],
                drug = unlist(strsplit(int_drug, "_"))[2]) %>%
  dplyr::select(c(time:resistance, efficacy:drug)) %>%
  pivot_wider(names_from = interval, values_from = efficacy)

drug_res$drug <- as.factor(drug_res$drug)
drug_res$resistance <- factor(drug_res$resistance, levels = c("Triple", "Quintuple"))

data_pt <- data.frame(time = c(28, 42, 28),
                      efficacy = c(0.94, 0.81, 0.90),
                      resistance = c("Triple", "Quintuple", "Quintuple"),
                      drug = rep("spaq", 3))
data_pt$resistance <- factor(data_pt$resistance, levels = c("Triple", "Quintuple"))

## pick up from here
ggplot(drug_res, aes(x = time, col = drug, fill = drug)) + 
  geom_line(aes(y = med), lwd = 0.8) + geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.2) + 
  facet_grid(resistance ~ .) +
  xlim(c(0, 60)) + 
  labs(x = "Time (days)", y = "Protective efficacy") + 
  geom_point(data = data_pt, aes(x = time, y = efficacy))
ggsave("output/figure_4.png", dpi = 300, width = 20, height = 15, units = "cm")


##------------------------------------------------------------------------------------------------------------
## create and save figures
# figure 1 -- step fit
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

ggsave("output/figure_1.png", dpi = 500, width = 20, height = 14, units = "cm")

# figure 2 - kaplan meier
control_q95 <- t(apply(control_mat, 2, quantile_95))
colnames(control_q95) <- sprintf("control_%s", colnames(control_q95))

treat_q95 <- t(apply(treat_mat, 2, quantile_95))
colnames(treat_q95) <- sprintf("treat_%s", colnames(treat_q95))

df_model <- data.frame(time = 0:(nrow(control_q95) - 1)) %>%
  bind_cols(control_q95) %>%
  bind_cols(treat_q95)

# --------------------------
# plot results

# plot against KM data
dat_trial %>%
  group_by(treat_arm) %>%
  dplyr::summarise(time = time.1,
                   n_infected = cumsum(n_infected)) %>%
  mutate(treat_arm = c("Control", "Treatment")[treat_arm],
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
  xlab("Time (days)") + ylab("Number infected") + 
  guides(col = guide_legend("Treatment arm"), fill = guide_legend("Treatment arm"))
ggsave("output/figure_2.png", dpi = 500, width = 20, height = 14, units = "cm")

## figure 3 - multipanel plot
pd_plot <- ggplot() +
  geom_line(data = pd_cri, aes(x = concentration, y = efficacy, group = sample),
            alpha = 0.1, col = "darkgrey") +
  geom_line(data = pd_median, aes(x = concentration, y = efficacy), lwd = 1) +
  scale_x_log10() + theme_bw() +
  labs(x = "log concentration (\u03bcg / ml)", y = "Protective efficacy")

pe_plot <- ggplot() + 
  geom_line(data = mcmc_post, aes(x = time/24, y = efficacy, group = sample),
            alpha = 0.05, col = "darkgrey") +
  theme_bw() +
  geom_line(data = mcmc_median, aes(x = time/24, y = efficacy), lwd = 1) +
  xlim(c(0.5/24, 60)) + labs(x = "Time (days)", y = "Protective efficacy")

pk_age <- ggplot() + 
  geom_line(data = pk_age_med, aes(x = time/24, y = median/1000, col = age_group)) + 
  theme_bw() + xlim(c(0,50)) + 
  geom_line(data = quadrature_pk_cri, aes(x = time/24, y = median/1000), 
            col = "black", linetype = 2) + 
  geom_ribbon(data = quadrature_pk_cri, aes(x = time/24, y = median/1000, 
                                            ymin = lower_cri/1000, ymax = upper_cri/1000),
              alpha = .1) +
  labs(x = "Time (days)", y = "Median sulfadoxine concentration (\u03bcg / ml)") +
  guides(col = guide_legend(title = "Age group")) + #, nrow=2, byrow=TRUE))  + 
  theme(legend.direction = "horizontal")

pk_nut <- ggplot() + 
  geom_line(data = pk_nut_med, aes(x = time/24, y = median/1000, col = z_score)) + 
  theme_bw() + xlim(c(0,50)) + 
  geom_line(data = quadrature_pk_cri, aes(x = time/24, y = median/1000), 
            col = "black", linetype = 2) + 
  geom_ribbon(data = quadrature_pk_cri, aes(x = time/24, y = median/1000, 
                                            ymin = lower_cri/1000, ymax = upper_cri/1000),
              alpha = .1) +
  labs(x = "Time (days)", y = "Median sulfadoxine concentration (\u03bcg / ml)") +
  guides(col = guide_legend(title = "Nutrition status")) + 
  theme(legend.direction = "horizontal")

pe_age <- ggplot() + geom_line(data = pe_sample, aes(x = time/24, y = efficacy, group = sample), 
                               alpha = 0.15, col = "darkgrey") +
  ylim(c(0,1)) + theme_bw() + xlim(c(0.5, 60)) +
  geom_line(data = pk_age_med, aes(x = time/24, y = efficacy, col = age_group), linewidth = 1) + 
  labs(x = "Time (days)", y = "Protective efficacy") +
  guides(col = guide_legend(title = "Age group")) 

pe_nut <- ggplot() + geom_line(data = pe_sample, aes(x = time/24, y = efficacy, group = sample), 
                     alpha = 0.15, col = "darkgrey") +
  ylim(c(0,1)) + theme_bw() + xlim(c(0.5, 60)) +
  geom_line(data = pk_nut_med, aes(x = time/24, y = efficacy, col = z_score), linewidth = 1) + 
  labs(x = "Time (days)", y = "Protective efficacy") +
  guides(col = guide_legend(title = "Nutrition status"))

get_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend) }

age_legend <- get_legend(pk_age)
nut_legend <- get_legend(pk_nut)

left_col <- ggpubr::ggarrange(pd_plot, pe_plot, ncol = 1, nrow = 2,
                              labels = c("A", "B"))
age_row <- ggpubr::ggarrange(pk_age + theme(legend.position = "none"), 
                             pe_age  + theme(legend.position = "none"),
                             ncol = 2,
                             labels = c("C", "E"))
age_row_legend <- ggpubr::ggarrange(age_row, age_legend, nrow = 2,
                                    heights = c(10, 1))
nut_row <- ggpubr::ggarrange(pk_nut + theme(legend.position = "none"), 
                             pe_nut  + theme(legend.position = "none"),
                             ncol = 2,
                             labels = c("D", "F"))
nut_row_legend <- ggpubr::ggarrange(nut_row, nut_legend, nrow = 2,
                                    heights = c(10, 1))

age_nut <- ggpubr::ggarrange(age_row_legend, nut_row_legend, nrow = 2)
fig_3 <- ggpubr::ggarrange(left_col, age_nut, ncol = 2, widths = c(1, 2))

jpeg("output/figure_3.jpg", width = 25, height = 15, units = "cm", res = 700)
ggpubr::ggarrange(left_col, age_nut, ncol = 2, widths = c(1, 2))
dev.off()
