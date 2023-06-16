## 04.posterior_estimates
mcmc <- readRDS("ignore/outputs/mcmc_raw_final.rds")

mcmc_post <- mcmc$output %>%
  dplyr::filter(phase == "sampling") %>% 
  tidyr::pivot_longer(cols = lambda_1:hill_power, names_to = "parameter", values_to = "value") %>% 
  dplyr::group_by(parameter)  %>%
  dplyr::mutate(scaling = if_else(grepl("lambda", parameter), 1000, 1)) %>%
  dplyr::summarise(scaling = scaling, 
                   median = median(value) * scaling,
                   lower_cri = quantile(value, probs = 0.025) * scaling,
                   upper_cri = quantile(value, probs = 0.975) * scaling) %>%
  distinct() %>%
  dplyr::arrange(factor(parameter, 
                        levels = c("min_prob", "half_point", "hill_power",
                                   paste0("lambda_", 1:13))))

## diagnostics
mcmc$diagnostics
