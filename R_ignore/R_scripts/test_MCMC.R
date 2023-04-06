
# read in drug concentration data
dat_raw <- readRDS("R_ignore/data/quadrature_pk_abd.rds")

# get weight of each group*quadrature combination
w_combined <- dat_raw %>%
  dplyr::filter(time == 0) %>%
  mutate(w = weighting * pop_prop) %>%
  pull(w)

# define parameters dataframe
df_params <- define_params(name = "FOI", min = 0, max = Inf,
                           name = "min_prob", min = 0, max = 1,
                           name = "half_point", min = 0, max = Inf,
                           name = "hill_power", min = 0, max = Inf)

# run MCMC
mcmc <- run_mcmc(data = list(ind = dat_raw$individual,
                             time = dat_raw$time,
                             drug_conc = dat_raw$drug_value.conc,
                             ind_weight = w_combined),
                 df_params = df_params,
                 burnin = 1e1,
                 samples = 1e2,
                 chains = 1)

plot_par(mcmc, show = "FOI")
