## Import libraries, datasets and functions
source("R/import_library.R")
source("R/import_covariates.R")
# Import list of simulations
list_all_sim <- readRDS("Output/simulation_set.RDS")
## Modify the interpretcontrol function from the surveillance package to 
## integrate the transmission potential to hhh4 models.
# In the new version stsobj@observed corresponds to the transmission potential,
# whereas control$response corresponds to the number of cases per region at 
# each time step. The transmission potential is used to computed the mean 
# number of expected cases at t, whereas the "response" data frame is used to
# generate the likelihood (i.e. it corresponds to the daily data).
assignInNamespace("interpretControl", my.interpretControl, ns = "surveillance", 
                  pos = "package:surveillance")

#### Apply report rate ####

report <- .7
list_all_sim_reported <- list_all_sim
set.seed(1)
list_all_sim_reported[[1]] <- lapply(list_all_sim[[1]], function(X){
  # Apply report rate
  data_i <- report_rate(data = X, report_rate = report)
  return(data_i)
})

#### Analysis daily data #### 

models_daily <- 
  function_hhh4_daily(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                      area_mat = area_ts, cov_mat = cov_ts, thresh = NA, 
                      distance_matrix = dist_mat_cent, prop_gen1 = .5, 
                      mean_si = 11.7, sd_si = 2.0, max_si = 50,
                      fun_wei = W_exp_gravity_tot)
# Save all daily models
saveRDS(models_daily, "Output/models_daily.RDS")

#### Analysis aggregated data ####

models_aggre <- 
  function_hhh4_aggreg(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                       area_mat = area_ts, cov_mat = cov_ts, thresh = NA, 
                       len_agg = 10, distance_matrix = dist_mat_cent, 
                       fun_wei = W_exp_gravity_tot)

# Save all aggregated models
saveRDS(models_aggre, "Output/models_aggregated.RDS")

