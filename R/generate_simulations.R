##### Import functions,libraries, and covariates

source("R/import_library.R")
source("R/import_covariates.R")

##### Generate simulations

list_params <- readRDS("Data/params_model_day_tot.RDS")

w_dens <- serial_interval(prop_gen1 = 1, mean_w1 = 11.7, sd_w1 = 2.0)
set.seed(1)
list_all_sim <- generate_all_sim(n_sim = 100, w_dens = w_dens, cov = cov_ts, 
                                 pop = pop_ts, distance_matrix = dist_mat_cent, 
                                 list_params = list_params, area = area_ts)

saveRDS(list_all_sim, "Output/simulation_set.RDS")
