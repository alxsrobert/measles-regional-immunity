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
map <- generate_map()

### Apply report rate

report <- .7
list_all_sim_reported <- list_all_sim
set.seed(1)
list_all_sim_reported[[1]] <- lapply(list_all_sim[[1]], function(X){
  # Apply report rate
  data_i <- report_rate(data = X, report_rate = report)
  return(data_i)
})

### Plots data
p <- plot_simulations(list_all_sim = list_all_sim, which_plots = c(45, 64), 
                      list_all_sim_reported = list_all_sim_reported, map = map)

### Analysis hhh4 runs 
# Generate figures, select which parameter to plot in variable "which"
plot_analysis(list_all_sim = list_all_sim_reported, hhh4_day = models_daily, 
              hhh4_agg = models_aggre, CI = .95, 
              which_plot = c("vacc_ar", "vacc_ne", "vacc_en",
                             "cat1_imm_ar", "cat1_imm_ne", "cat1_imm_en",
                             "cat2_imm_ar", "cat2_imm_ne", "cat2_imm_en"), 
              mains = c("Autoregressive - Vaccination", "Neighbourhood - Vaccination",
                        "Endemic - Vaccination", "Autoregressive - Incidence 1", 
                        "Neighbourhood - Incidence 1", "Endemic - Incidence 1",
                        "Autoregressive - Incidence 2", "Neighbourhood - Incidence 2",
                        "Endemic - Incidence 2"),
              exclude_overdisp = T)

### Proportion in each component
plot_prop_comp(list_all_sim = list_all_sim, hhh4_day = models_daily, 
               hhh4_agg = models_aggre)
