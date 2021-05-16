## Import libraries, datasets and functions
source("R/import_library.R")
source("R/import_covariates.R")

library(tictoc)

# Set the interpretcontrol function
assignInNamespace("interpretControl", my.interpretControl, ns = "surveillance", 
                  pos = "package:surveillance")
models_daily <- readRDS("Output/models_daily.RDS")
models_aggre <- readRDS("Output/models_aggregated.RDS")
dates_aggreg <- seq(1, nrow(models_daily[[1]]$control$data$response), 10)
all_dates <- dates_aggreg[length(dates_aggreg) - 35:1] - 1
all_dates_agg <- round((dates_aggreg[length(dates_aggreg) - 35:1])/10)

##### Generate the calibration scores  

scores_day <- calibration_model(model_list = models_daily[1:25], daily = T, 
                                all_dates = all_dates, w_dens = w_dens, 
                                new_fit = T, period_pred = 10, k = 20000)
scores_agg <- calibration_model(model_list = models_aggre, daily = F, 
                                all_dates = all_dates_agg, new_fit = T, 
                                w_dens = w_dens)

saveRDS(scores_day, "Output/scores_day3.RDS")
saveRDS(scores_agg, "Output/scores_aggreg.RDS")
