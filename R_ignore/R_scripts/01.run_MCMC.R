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
dat_drug <- readRDS("R_ignore/data/quadrature_pk.rds")

# get weight of each group*quadrature combination
w_combined <- dat_drug %>%
  dplyr::filter(time == 0) %>%
  mutate(w = weighting * pop_prop) %>%
  pull(w)
w_combined <- w_combined / sum(w_combined)

# read in trial data and split into control vs. treatment
dat_trial <- read.csv("R_ignore/data/cisse_data.csv")
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
                             data_treat = dat_treat),
                 burnin = 1e2,
                 samples = 1e2,
                 chains = 5)


# --------------------------
# exploratory plots

mcmc$output %>%
  #dplyr::filter(phase == "sampling") %>%
  dplyr::select(-c(phase, logprior, loglikelihood)) %>%
  pivot_longer(cols = -c(chain, iteration), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = factor(parameter, levels = c(sprintf("lambda_%s", 1:13), "min_prob", "half_point", "hill_power"))) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = iteration, y = value, col = chain), size = 0.5) +
  facet_wrap(~parameter, scales = "free_y")

# param_bad <- mcmc$output %>%
#   dplyr::filter(chain == 3) %>%
#   dplyr::filter(iteration == 1100) %>%
#   dplyr::select(-c(chain, phase, iteration, logprior, loglikelihood)) %>%
#   unlist()

get_loglike(data = list(data_drug = dat_drug,
                        ind_weight = w_combined,
                        data_control = dat_control,
                        data_treat = dat_treat),
            params = as.list(param_bad))

param_prop <- param_bad
param_prop["half_point"] <- 0.1

get_loglike(data = list(data_drug = dat_drug,
                        ind_weight = w_combined,
                        data_control = dat_control,
                        data_treat = dat_treat),
            params = as.list(param_prop))

# --------------------------
# save output to file

if (FALSE) {
  saveRDS(mcmc, file = "R_ignore/outputs/mcmc_raw.rds")
}


