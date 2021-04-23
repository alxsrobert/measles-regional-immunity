## Import libraries, datasets and functions
source("R/import_library.R")
source("R/import_covariates.R")
# Import list of simulations
list_all_sim <- readRDS("Output/simulation_set.RDS")
# Set the number interpretcontrol function
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
                      area_mat = area_ts, cov_mat = cov_ts, thresh = 45, 
                      distance_matrix = dist_mat_cent, prop_gen1 = .5, 
                      mean_si = 11.7, sd_si = 2.0, max_si = 50,
                      report_rate = .7)
# Save all daily models
saveRDS(models_daily, "Output/models_daily.RDS")

#### Analysis aggregated data ####

models_aggre <- 
  function_hhh4_aggreg(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                       area_mat = area_ts, cov_mat = cov_ts, thresh = 45, 
                       len_agg = 10, distance_matrix = dist_mat_cent,
                       report_rate = .7)

# Save all aggregated models
saveRDS(models_aggre, "Output/models_aggregated.RDS")

