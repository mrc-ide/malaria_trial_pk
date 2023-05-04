# 01.run_MCMC.R
#
# Author: Bob Verity
# Date: 2023-04-14
#
# Purpose:
# Reads in drug concentration data from PK model and trial data and runs MCMC to
# estimate epi and PD model parameters. Saves MCMC output to file for
# summarising elsewhere.
#
# ------------------------------------------------------------------

# load packages
library(tidyverse)

# load this package
devtools::load_all(".")

# read in drug concentration data
dat_drug <- readRDS("data/quadrature_pk.rds")

# get weight of each group*quadrature combination (same for all time points)
w_combined <- dat_drug %>%
  dplyr::filter(time == 0) %>%
  mutate(w = weighting * pop_prop) %>%
  pull(w)

# get EIR adjustment (same for all time points)
eir_adjustment <- dat_drug %>%
  dplyr::filter(time == 0) %>%
  pull(eir_adjustment)

# read in trial data and split into control vs. treatment
dat_trial <- read.csv("data/cisse_data.csv")
dat_control <- dat_trial %>%
  dplyr::filter(treat_arm == 1) %>%
  dplyr::select(n_patients, n_infected, time, time.1)
dat_treat <- dat_trial %>%
  dplyr::filter(treat_arm == 2) %>%
  dplyr::select(n_patients, n_infected, time, time.1)

# --------------------------

# run MCMC
set.seed(1)
mcmc <- run_mcmc(data = list(data_drug = dat_drug,
                             ind_weight = w_combined,
                             data_control = dat_control,
                             data_treat = dat_treat,
                             eir_adjustment = eir_adjustment),
                 burnin = 1e3,#1e3,
                 samples = 1e3,#5e3,
                 chains = 1)#10)


# --------------------------
# exploratory plots

# trace plots
mcmc$output %>%
  dplyr::filter(phase == "sampling") %>%
  dplyr::select(-c(phase, logprior, loglikelihood)) %>%
  pivot_longer(cols = -c(chain, iteration), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c(sprintf("lambda_%s", 1:13), "min_prob", "half_point", "hill_power"))) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value, col = chain), size = 0.5) +
  facet_wrap(~parameter, scales = "free_y")

# calculate quantiles
df_quantile <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>%
  dplyr::select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
  apply(2, quantile_95) %>%
  t() %>%
  as.data.frame()

# quantile plots
df_quantile$names <- row.names(df_quantile)
df_quantile$is_lambda <- grepl("lambda", df_quantile$names)
plot1 <- df_quantile %>%
  dplyr::filter(is_lambda == 1) %>%
  ggplot() + theme_bw() +
  geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5, y = names)) +
  xlab("Value") + xlab("Parameter")
plot2 <- df_quantile %>%
  dplyr::filter(is_lambda == 0) %>%
  ggplot() + theme_bw() +
  geom_errorbar(aes(xmin = Q2.5, xmax = Q97.5, y = names)) +
  xlab("Value") + xlab("Parameter")
cowplot::plot_grid(plot1, plot2)

# --------------------------
# save output to file

if (FALSE) {
  saveRDS(mcmc, file = "ignore/outputs/mcmc_raw.rds")
}


