## Import libraries, datasets and functions
source("R/import_library.R")
source("R/import_covariates.R")

## Modify the interpretcontrol function from the surveillance package to 
## integrate the transmission potential to hhh4 models.
# In the new version stsobj@observed corresponds to the transmission potential,
# whereas control$response corresponds to the number of cases per region at 
# each time step. The transmission potential is used to computed the mean 
# number of expected cases at t, whereas the "response" data frame is used to
# generate the likelihood (i.e. it corresponds to the daily data).
assignInNamespace("interpretControl", my.interpretControl, ns = "surveillance", 
                  pos = "package:surveillance")
# Import daily models
models_daily <- c(readRDS("Output/models_daily1.RDS"),
                  readRDS("Output/models_daily2.RDS"),
                  readRDS("Output/models_daily3.RDS"),
                  readRDS("Output/models_daily4.RDS"))
# import aggregated models
models_aggre <- readRDS("Output/models_aggregated.RDS")
# Generate dates of aggregation
dates_aggreg <- seq(1, nrow(models_daily[[1]]$control$data$response), 10)
all_dates <- dates_aggreg[length(dates_aggreg) - 35:1] - 1
all_dates_agg <- round((dates_aggreg[length(dates_aggreg) - 35:1])/10)

##### Generate the calibration scores  

scores_day <- calibration_model(model_list = models_daily, daily = T, 
                                all_dates = all_dates, w_dens = w_dens, 
                                new_fit = T, period_pred = 10, k = 20000)
scores_agg <- calibration_model(model_list = models_aggre, daily = F, 
                                all_dates = all_dates_agg, new_fit = T, 
                                w_dens = w_dens)

saveRDS(scores_day, "Output/scores_day.RDS")
saveRDS(scores_agg, "Output/scores_aggreg.RDS")
