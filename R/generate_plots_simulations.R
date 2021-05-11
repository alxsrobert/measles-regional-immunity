## Import libraries, datasets and functions
source("R/import_library.R")
source("R/import_covariates.R")
# Import simulations
list_all_sim <- readRDS("Output/simulation_set.RDS")
# Import model fits (daily and aggregated)
models_daily <- c(readRDS("Output/models_daily1.RDS"),
                  readRDS("Output/models_daily2.RDS"),
                  readRDS("Output/models_daily3.RDS"),
                  readRDS("Output/models_daily4.RDS"))
models_aggre <- readRDS("Output/models_aggregated.RDS")

### Apply report rate

report <- .7
list_all_sim_reported <- list_all_sim
set.seed(1)
list_all_sim_reported[[1]] <- lapply(list_all_sim[[1]], function(X){
  # Apply report rate
  data_i <- report_rate(data = X, report_rate = report)
  return(data_i)
})


### Analysis hhh4 runs 
# Generate figures, select which parameter to plot in variable "which"
plot_analysis(list_all_sim = list_all_sim_reported, hhh4_day = models_daily, 
              hhh4_agg = models_aggre, CI = .95, 
              which_plot = c("vacc_ar", "vacc_ne", "vacc_en",
                             "cat1_imm_ar", "cat1_imm_ne", "cat1_imm_en",
                             "cat2_imm_ar", "cat2_imm_ne", "cat2_imm_en"), 
              exclude_overdisp = T)

