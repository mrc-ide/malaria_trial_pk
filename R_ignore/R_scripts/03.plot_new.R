## plotting script from scratch to try and fix the issues with the current plotting script
library(tidyverse)
library(pammtools)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra) 
library(reldist)
library(ggplot2)

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

source("R_ignore/R_scripts/plot_functions.R")

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
n_sub <- 1e3
mcmc_sub <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>%
  sample_n(n_sub)

# figure 1
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

arm_labs <- c("Control", "SP+AS")
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

rate_plot <- ggplot() + geom_step(data = df_final, 
                                  aes(x = time/24, y = rate * 1000, 
                                      col = treat_arm)) + 
  facet_grid(facet_factor ~ ., scales = "free_y") +
  theme_bw() + theme(strip.text.y = element_blank()) +
  geom_stepribbon(data = mcmc_summary_final, 
                  aes(x = time/24, ymin = rate_2.5* 1000, ymax = rate_97.5* 1000, fill = treat_arm), alpha = 0.2) + 
  geom_step(data = mcmc_summary_final, aes(x = time/24, y = rate_50* 1000, col = treat_arm), linetype = 2) + 
  scale_color_discrete(labels = arm_labs) + 
  scale_fill_discrete(labels = arm_labs) +
  guides(fill = guide_legend(title = "Treatment arm"),
         col = guide_legend(title = "Treatment arm")) +
  labs(x = "Time (days)", y = "Infection rate per 1,000 \nchildren at risk per day ") + 
  geom_vline(xintercept = c(0, 28, 56), linetype = "dotted")

# adapt so that the plotting colours and labels are consistent
df_model <- readRDS("data/df_model.rds")
df_model_control <- df_model %>% dplyr::select(time:control_Q97.5)
df_model_spas <- df_model %>% dplyr::select(time, treat_Q2.5:treat_Q97.5)
df_model_control$arm <- "Control"
df_model_spas$arm <- "SP+AS"
names(df_model_control) <- c("time", "Q2.5", "Q50", "Q97.5", "arm")
names(df_model_spas) <- c("time", "Q2.5", "Q50", "Q97.5", "arm")
df_model_long <- rbind(df_model_control, df_model_spas)
df_model_long$arm <- factor(df_model_long$arm)

km_plot <- dat_trial %>%
  group_by(treat_arm) %>%
  dplyr::summarise(time = time.1,
                   n_infected = cumsum(n_infected)) %>%
  mutate(treat_arm = c("Control", "SP+AS")[treat_arm],
         treat_arm = factor(treat_arm)) %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = time / 24, y = Q50, col = arm), 
            alpha = 0.2, data = df_model_long, lty = 2, lwd = 0.8) +
  geom_ribbon(aes(x = time / 24, ymin = Q2.5,
                  ymax = Q97.5,fill = arm,),
              alpha = 0.2, data = df_model_long) +
  #geom_step(aes(x = time / 24, y = n_infected, col = treat_arm), direction = "hv") +
  geom_point(aes(x = time / 24, y = n_infected, col = treat_arm), size = 0.8) +
  xlab("Time (days)") + ylab("\nNumber infected") + 
  guides(col = guide_legend("Treatment arm"), fill = guide_legend("Treatment arm")) + 
  geom_vline(xintercept = c(0, 28, 56), linetype = "dotted")

ggpubr::ggarrange(km_plot, rate_plot, ncol = 1, nrow = 2,
                  labels = c("A", "B"), common.legend = TRUE)
ggsave("output/figure_1.png", dpi = 500, width = 20, height = 20, units = "cm")

### figure 2
# A shows the median PD and posterior samples
# B shows the median PE and posterior samples
# C shows the median PK and credible interval + median by age
# D shows the median PK and credible interval + median by nutrition
# E shows the mean PE and posterior samples + mean PE by age
# F shows the mean PE and posterior samples + mean PE by nutrition
## pharamacodynamics
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
concentration_vec <- seq(0.1, 150, by = 0.1)
pd_median <- median_pd(conc_vec = concentration_vec,
                       mcmc_cri = mcmc_cri)
pd_cri <- sample_pd(mcmc_output = mcmc$output, num_samples = 1000,
                    conc_vec = concentration_vec)

plot_2a <- ggplot() +
  geom_line(data = pd_cri, aes(x = concentration, y = efficacy, group = sample),
            alpha = 0.1, col = "darkgrey") +
  geom_line(data = pd_median, aes(x = concentration, y = efficacy), lwd = 1) +
  scale_x_log10() + theme_bw() +
  labs(x = "log median concentration (\u03bcg / ml)", y = "Median protective efficacy")

# median PE for now
median_df <- median_pk(quadrature_pk_1)
mcmc_post <- sample_posterior(mcmc_output = filter(mcmc$output, phase == "sampling"), 
                              num_samples = 1000, 
                              median_df = median_df)
head(mcmc_post)

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

mcmc_median <- median_posterior(mcmc_summary = mcmc_cri, median_df = median_df)

plot_2b <- ggplot() + 
  geom_line(data = mcmc_post, aes(x = time/24, y = efficacy, group = sample),
            alpha = 0.05, col = "darkgrey") +
  theme_bw() +
  geom_line(data = mcmc_median, aes(x = time/24, y = efficacy), lwd = 1) +
  xlim(c(0.5/24, 60)) + labs(x = "Time (days)", y = "Median concentration (\u03bcg / ml)") + 
  geom_vline(xintercept = 35, lty = 2, lwd = 0.4)

# pk plots
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
                                            weight = weighting * pop_prop, q = 0.975),
                   mean = wtd.mean(drug_value.conc,
                                   weight = weighting * pop_prop))
pk_nut_med <- quadrature_pk_age_nut %>% 
  dplyr::group_by(time, z_score) %>% 
  dplyr::summarise(lower_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.025),
                   median = median(drug_value.conc),
                   upper_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.975), 
                   mean = wtd.mean(drug_value.conc,
                                   weight = weighting * pop_prop))

quadrature_pk_cri <- quadrature_pk_1 %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(lower_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.025),
                   median = median(drug_value.conc),
                   upper_cri = wtd.quantile(drug_value.conc, 
                                            weight = weighting * pop_prop, q = 0.975),
                   mean = wtd.mean(drug_value.conc,
                                   weight = weighting * pop_prop))

pk_age_med$age_group <- factor(pk_age_med$age_group)
pk_nut_med$z_score <- factor(pk_nut_med$z_score)

plot_2c <- ggplot() + 
  geom_line(data = pk_age_med, aes(x = time/24, y = median, col = age_group)) + 
  theme_bw() + xlim(c(0,50)) + theme(legend.direction = "horizontal", legend.position = "bottom") + 
  geom_line(data = quadrature_pk_cri, aes(x = time/24, y = median), 
            col = "black", linetype = 2) + 
  geom_ribbon(data = quadrature_pk_cri, aes(x = time/24, y = median, 
                                            ymin = lower_cri, ymax = upper_cri),
              alpha = .1) +
  labs(x = "days", y = "Median sulfadoxine \nconcentration (\u03bcg / ml)") +
  guides(col = guide_legend(title = "Age group"))

plot_2d <- ggplot() + 
  geom_line(data = pk_nut_med, aes(x = time/24, y = median, col = z_score)) + 
  theme_bw() + xlim(c(0,50)) + theme(legend.direction = "horizontal", legend.position = "bottom") + 
  geom_line(data = quadrature_pk_cri, aes(x = time/24, y = median), 
            col = "black", linetype = 2) + 
  geom_ribbon(data = quadrature_pk_cri, aes(x = time/24, y = median, 
                                            ymin = lower_cri, ymax = upper_cri),
              alpha = .1) +
  labs(x = "days", y = "Median sulfadoxine \nconcentration (\u03bcg / ml)") +
  guides(col = guide_legend(title = "Nutritional status"))

## mean PE by status
mcmc_output <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>%
  dplyr::select(c(min_prob, half_point, hill_power)) %>%
  sample_n(size = 50) # update to 250 for final run

# pe_sample uses weighted means i.e. gives weighted mean per iteration
t0 <- Sys.time()
pe_sample <- pe_cri(mcmc_output, quadrature_pk_1)
t0 - Sys.time()

# make subpops by age and nutrition
head(pk_age_med)
head(pk_nut_med)

for(i in unique(pk_age_med$age_group)) {
  assign(paste0("pk_age_", gsub(" years", "", x = gsub(" - ", "_", x = i))), 
         subset(quadrature_pk_age_nut, age_group == i))
}

for(i in unique(pk_nut_med$z_score)) {
  assign(substr(paste0("pk_nut_", x = gsub(" - ", "_", x = i)), 1, 10), 
         subset(quadrature_pk_age_nut, z_score == i))
}

t0 <- Sys.time()
pk_age_0_1_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_age_0_1) 
Sys.time() - t0


pk_age_1_2_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_age_1_2) 
pk_age_2_3_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_age_2_3) 
pk_age_3_4_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_age_3_4) 
pk_age_4_5_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_age_4_5) 

pk_nut_mal_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_nut_mal) 
pk_nut_not_sample <- pe_cri(mcmc_output = mcmc_output, quad_1 = pk_nut_not)

pk_age_0_1_sample$age_group <- "0 - 1 years"
pk_age_1_2_sample$age_group <- "1 - 2 years"
pk_age_2_3_sample$age_group <- "2 - 3 years"
pk_age_3_4_sample$age_group <- "3 - 4 years"
pk_age_4_5_sample$age_group <- "4 - 5 years"
pk_nut_mal_sample$z_score <- "malnourished"
pk_nut_not_sample$z_score <- "not malnourished"

pk_age_0_1_sample <- pk_age_0_1_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))
pk_age_1_2_sample <- pk_age_1_2_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))
pk_age_2_3_sample <- pk_age_2_3_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))
pk_age_3_4_sample <- pk_age_3_4_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))
pk_age_4_5_sample <- pk_age_4_5_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))
pk_nut_mal_sample <- pk_nut_mal_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))
pk_nut_not_sample <- pk_nut_not_sample %>% dplyr::group_by(time) %>%
  dplyr::mutate(mean = mean(efficacy))

#save these
saveRDS(pk_age_0_1_sample, "data/pk_age_0_1_sample.rds")
saveRDS(pk_age_1_2_sample, "data/pk_age_1_2_sample.rds")
saveRDS(pk_age_2_3_sample, "data/pk_age_2_3_sample.rds")
saveRDS(pk_age_3_4_sample, "data/pk_age_3_4_sample.rds")
saveRDS(pk_age_4_5_sample, "data/pk_age_4_5_sample.rds")
saveRDS(pk_nut_mal_sample, "data/pk_nut_mal_sample.rds")
saveRDS(pk_nut_not_sample, "data/pk_nut_not_sample.rds")

# combine into one
pk_age_mean <- rbind(pk_age_0_1_sample,
                     pk_age_1_2_sample,
                     pk_age_2_3_sample,
                     pk_age_3_4_sample,
                     pk_age_4_5_sample)

pk_nut_mean <- rbind(pk_nut_mal_sample, 
                     pk_nut_not_sample)

pk_age_mean <- pk_age_mean %>% dplyr::group_by(time) %>% 
  dplyr::mutate(pop_mean = mean(efficacy)) #%>%
# dplyr::filter(iteration <= 10)
pk_nut_mean <- pk_nut_mean %>% dplyr::group_by(time) %>% 
  dplyr::mutate(pop_mean = mean(efficacy)) #%>%
# dplyr::filter(iteration <= 25)
pk_age_mean$age_group <- pk_age_mean$age_group
pk_nut_mean$z_score <- pk_nut_mean$z_score

pk_age_mean$unique <- paste0(pk_age_mean$iteration, "_", pk_age_mean$age_group)
pk_nut_mean$unique <- paste0(pk_nut_mean$iteration, "_", pk_nut_mean$z_score)

plot_2e <- ggplot(pk_age_mean) + 
  geom_line(aes(x = time/24, y = efficacy, group = unique, col = age_group), alpha = 0.05) +
  geom_line(aes(x = time/24, y = pop_mean), lwd = 1) +
  geom_line(aes(x = time/24, y = mean, col = age_group), lwd = 1) +
  ylim(c(0,1)) + theme_bw() + xlim(c(0.5, 60)) +
  guides(col = guide_legend(title = "Age group")) +
  labs(x = "Time (days)", y = "Mean protective efficacy") + 
  geom_vline(xintercept = 35, lty = 2, lwd = 0.4)

plot_2f <- ggplot(pk_nut_mean) + 
  geom_line(aes(x = time/24, y = efficacy, group = unique, col = z_score), alpha = 0.05) +
  geom_line(aes(x = time/24, y = pop_mean), lwd = 1) +
  geom_line(aes(x = time/24, y = mean, col = z_score), lwd = 1) +
  ylim(c(0,1)) + theme_bw() + xlim(c(0.5, 60)) +
  guides(col = guide_legend(title = "Nutritional status")) +
  labs(x = "Time (days)", y = "Mean protective efficacy") + 
  geom_vline(xintercept = 35, lty = 2, lwd = 0.4)

get_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend) }

age_legend <- get_legend(plot_2c)
nut_legend <- get_legend(plot_2d)

left_col <- ggpubr::ggarrange(plot_2a, plot_2b, ncol = 1, nrow = 2,
                              labels = c("A", "B"))
age_row <- ggpubr::ggarrange(plot_2c + theme(legend.position = "none"), 
                             plot_2e  + theme(legend.position = "none"),
                             ncol = 2,
                             labels = c("C", "E"))
age_row_legend <- ggpubr::ggarrange(age_row, age_legend, nrow = 2,
                                    heights = c(10, 1))
nut_row <- ggpubr::ggarrange(plot_2d + theme(legend.position = "none"), 
                             plot_2f  + theme(legend.position = "none"),
                             ncol = 2,
                             labels = c("D", "F"))
nut_row_legend <- ggpubr::ggarrange(nut_row, nut_legend, nrow = 2,
                                    heights = c(10, 1))

age_nut <- ggpubr::ggarrange(age_row_legend, nut_row_legend, nrow = 2)
fig_2 <- ggpubr::ggarrange(left_col, age_nut, ncol = 2, widths = c(1, 2))

jpeg("output/figure_2.jpg", width = 30, height = 15, units = "cm", res = 700)
fig_2
dev.off()

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
drug_res <- drug_res %>%
  dplyr::mutate(drug = toupper(drug),
                drug = if_else(drug == "SPAQ", "SP+AQ", drug))

data_pt <- data.frame(time = c(28, 42, 28),
                      efficacy = c(0.94, 0.81, 0.90),
                      resistance = c("Triple", "Quintuple", "Quintuple"),
                      drug = rep("SP+AQ", 3))
data_pt$resistance <- factor(data_pt$resistance, levels = c("Triple", "Quintuple"))

ggplot(drug_res, aes(x = time, fill = drug)) + 
  geom_line(aes(y = med, col = drug), lwd = 0.8) + geom_ribbon(aes(ymin = low, ymax = upp), alpha = 0.2) + 
  facet_grid(resistance ~ .) +
  xlim(c(0, 60)) + 
  labs(x = "Time (days)", y = "Protective efficacy") + 
  geom_point(data = data_pt, aes(x = time, y = efficacy)) +
  guides(col = guide_legend("Drug"), fill = guide_legend("Drug"))
ggsave("output/figure_3.png", dpi = 300, width = 20, height = 15, units = "cm")

drug_res %>% dplyr::filter(resistance == "Quintuple") %>% 
  dplyr::filter(drug == "spaq") %>%
  dplyr::filter(time %in% seq(0, 28)) %>%
  dplyr::reframe(mean = mean(med),
                 lower_cri = mean(low),
                 upper_cri = mean(upp))
