#######################################################
######### Generate analysis on simulated data #########
#######################################################

## Import libraries, datasets and functions
source("R/import_library.R")
source("R/import_covariates.R")
source("R/function_Main_figures.R")
source("R/function_Supplement_figures.R")
source("R/function_sensitivity.R")
## Modify the interpretcontrol function from the surveillance package to 
## integrate the transmission potential to hhh4 models.
# In the new version stsobj@observed corresponds to the transmission potential,
# whereas control$response corresponds to the number of cases per region at 
# each time step. The transmission potential is used to computed the mean 
# number of expected cases at t, whereas the "response" data frame is used to
# generate the likelihood (i.e. it corresponds to the daily data).
assignInNamespace("interpretControl", my.interpretControl, ns = "surveillance", 
                  pos = "package:surveillance")
list_params <- readRDS("Data/params_model_day_tot.RDS")
w_comp <- serial_interval(prop_gen1 = .5, mean_w1 = 11.7, sd_w1 = 2.0)

#### Generate the simulated outbreaks ####

# Simulate 1 outbreak
set.seed(2)
list_all_sim <- generate_all_sim(n_sim = 1, w_dens = w_dens, cov = cov_ts, 
                                 pop = pop_ts, distance_matrix = dist_mat_cent, 
                                 list_params = list_params, area = area_ts, 
                                 rd_params = F)

# Apply report rate
report <- .7
list_all_sim_reported <- list_all_sim
set.seed(2)
list_all_sim_reported[[1]] <- lapply(list_all_sim[[1]], function(X){
  # Apply report rate
  data_i <- report_rate(data = X, report_rate = report)
  return(data_i)
})

#### Fit daily model using Neighbourhood and Total distance matrix ####

model_day_tot <- 
  function_hhh4_daily(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                      area_mat = area_ts, cov_mat = cov_ts, thresh = c(10, NA), 
                      distance_matrix = dist_mat_cent, prop_gen1 = .5, 
                      mean_si = 11.7, sd_si = 2.0, max_si = 50,
                      fun_wei = W_exp_gravity_tot)

model_day_nei <- 
  function_hhh4_daily(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                      area_mat = area_ts, cov_mat = cov_ts, thresh = c(10, NA), 
                      distance_matrix = dist_mat_nei, prop_gen1 = .5, 
                      mean_si = 11.7, sd_si = 2.0, max_si = 50,
                      fun_wei = W_exp_gravity_nei)

#### Simulate 1y predictions ####

days <- seq(as.Date("2017-12-01"), as.Date("2018-12-31"), 1)

sim_1y <- simulation_main_figures(model = model_day_tot[[1]], n_param = 10, 
                                  days = days, n_sim_param = 10, pop_mat = pop_ts,
                                  data = model_day_tot[[1]]$control$data$response,
                                  w_dens = w_comp)
sim_1y_nei <- simulation_main_figures(model = model_day_nei[[1]], n_param = 10, 
                                      days = days, n_sim_param = 10, pop_mat = pop_ts,
                                      data = model_day_tot[[1]]$control$data$response,
                                      w_dens = w_comp)

regions <- c("FR623", "FR624", "FR824", "FR101")

sim_loc <- list()
sim_loc_better <- list()
sim_loc_worse <- list()

sim_loc_nei <- list()

n_import <- 10
for(i in seq_along(regions)){
  reg_i <- regions[i]
  set.seed(1)
  dates_imp <- sample(x = days[days < "2018-01-01"], 
                      size = n_import, replace = T)
  
  sim_i <- 
    simulation_main_figures(model = model_day_tot[[1]], n_param = 10, 
                            days = days, n_sim_param = 10, pop_mat = pop_ts,
                            data = model_day_tot[[1]]$control$data$response,
                            w_dens = w_comp, reg_import = reg_i,
                            dates_import = dates_imp)
  sim_loc[[i]] <- sim_i[[1]]
  sim_loc_better[[i]] <- sim_i[[2]]
  sim_loc_worse[[i]] <- sim_i[[3]]
  
  sim_i_nei <- 
    simulation_main_figures(model = model_day_nei[[1]], n_param = 10, 
                            days = days, n_sim_param = 10, pop_mat = pop_ts,
                            data = model_day_nei[[1]]$control$data$response,
                            w_dens = w_comp, reg_import = reg_i, 
                            dates_import = dates_imp)
  sim_loc_nei[[i]] <- sim_i_nei[[1]]
}


#### Impact of categories of incidence ####

sens_incid_tot <- 
  sens_incid(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
             area_mat = area_ts, fun_wei = W_exp_gravity_tot, prop_gen1 = .5,
             distance_matrix = dist_mat_cent, mean_si = 11.7, sd_si = 2.0, 
             max_si = 50, cov_mat = cov_ts)
sens_incid_nei <- 
  sens_incid(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
             area_mat = area_ts, fun_wei = W_exp_gravity_nei, prop_gen1 = .5,
             distance_matrix = dist_mat_nei, mean_si = 11.7, sd_si = 2.0, 
             max_si = 50, cov_mat = cov_ts)

#### Sensitivity analysis serial interval ####

sens_si_tot <- sens_SI(vals = seq(0, 1, .1), start_coef = fixef(model_day_tot[[1]]),
                       list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                       area_mat = area_ts, cov_mat = cov_ts, thresh = c(10, NA), 
                       distance_matrix = dist_mat_cent, mean_si = 11.7,
                       sd_si = 2.0, max_si = 50, fun_wei = W_exp_gravity_tot)

sens_si_nei <- sens_SI(vals = seq(0, 1, .1), start_coef = fixef(model_day_nei[[1]]),
                       list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                       area_mat = area_ts, cov_mat = cov_ts, thresh = c(10, NA), 
                       distance_matrix = dist_mat_nei, mean_si = 11.7,
                       sd_si = 2.0, max_si = 50, fun_wei = W_exp_gravity_nei)

#### Sensitivity analysis vaccines ####

sens_vax_tot <- sens_vax(n = 10, start_coef = fixef(model_day_tot[[1]]),
                         list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                         area_mat = area_ts, thresh = c(10, NA), prop_gen1 = .5, 
                         distance_matrix = dist_mat_cent, mean_si = 11.7,
                         sd_si = 2.0, max_si = 50, fun_wei = W_exp_gravity_tot, 
                         corres = corres)

sens_vax_nei <- sens_vax(n = 10, start_coef = fixef(model_day_nei[[1]]),
                         list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                         area_mat = area_ts, thresh = c(10, NA), prop_gen1 = .5,
                         distance_matrix = dist_mat_nei, mean_si = 11.7,
                         sd_si = 2.0, max_si = 50, fun_wei = W_exp_gravity_nei, 
                         corres = corres)

#### Weekday effect ####

model_tot_weekday <- sens_weekday(list_all_sim = list_all_sim_reported, 
                                  pop_mat = pop_ts, area_mat = area_ts, 
                                  cov_mat = cov_ts, thresh = c(10, NA), prop_gen1 = .5,
                                  mean_si = 11.7, sd_si = 2.0, max_si = 50, 
                                  distance_matrix = dist_mat_cent, 
                                  fun_wei = W_exp_gravity_tot)


#### Calibration study ####

w_calib <- serial_interval(prop_gen1 = .5, mean_w1 = 11.7, sd_w1 = 2, 
                           max_days = 50)
dates_data <- as.Date(rownames(list_all_sim[[1]][[1]]))

dates_calib <- which(dates_data > as.Date("2016-01-03") &
                       dates_data < as.Date("2017-12-14"))
dates_calib_agg <- unique(trunc(dates_calib/10))
# 3 days ahead, Model 1
calib_3_tot <- calibration_model(model_list = model_day_tot, daily = T, 
                                 w_dens = w_calib, all_dates = dates_calib, 
                                 period_pred = 3)
# 7 days ahead, Model 1
calib_7_tot <- calibration_model(model_list = model_day_tot, daily = T, 
                                 w_dens = w_calib, all_dates = dates_calib, 
                                 period_pred = 7)
# 10 days ahead, Model 1
calib_10_tot <- calibration_model(model_list = model_day_tot, daily = T, 
                                  w_dens = w_calib, all_dates = dates_calib, 
                                  period_pred = 10)
# 14 days ahead, Model 1
calib_14_tot <- calibration_model(model_list = model_day_tot, daily = T, 
                                  w_dens = w_calib, all_dates = dates_calib, 
                                  period_pred = 14)
list_calib_tot <- list(calib_3_tot, calib_7_tot, calib_10_tot, calib_14_tot)

# 3 days ahead, Model 2
calib_3_nei <- calibration_model(model_list = model_day_nei, daily = T, 
                                 w_dens = w_calib, all_dates = dates_calib, 
                                 period_pred = 3)
# 7 days ahead, Model 2
calib_7_nei <- calibration_model(model_list = model_day_nei, daily = T, 
                                 w_dens = w_calib, all_dates = dates_calib, 
                                 period_pred = 7)
# 10 days ahead, Model 2
calib_10_nei <- calibration_model(model_list = model_day_nei, daily = T, 
                                  w_dens = w_calib, all_dates = dates_calib, 
                                  period_pred = 10)
# 14 days ahead, Model 2
calib_14_nei <- calibration_model(model_list = model_day_nei, daily = T, 
                                  w_dens = w_calib, all_dates = dates_calib, 
                                  period_pred = 14)
list_calib_nei <- list(calib_3_nei, calib_7_nei, calib_10_nei, calib_14_nei)

## Aggregated calibration 

models_aggre <- 
  function_hhh4_aggreg(list_all_sim = list_all_sim_reported, pop_mat = pop_ts, 
                       area_mat = area_ts, cov_mat = cov_ts, thresh = c(10, NA), 
                       len_agg = 10, distance_matrix = dist_mat_cent, 
                       fun_wei = W_exp_gravity_tot)
scores_agg <- calibration_model(model_list = models_aggre, daily = F, 
                                all_dates = dates_calib_agg, new_fit = T, 
                                w_dens = w_dens)



#### Save all files ####

# Group all important elements in a list
output_analysis <- list(data_rep = list_all_sim_reported[[1]][[1]], 
                        model_tot = model_day_tot[[1]],
                        model_nei = model_day_nei[[1]],
                        model_tot_weekday = model_tot_weekday,
                        sim_1y_tot = sim_1y,
                        sim_1y_nei = sim_1y_nei,
                        days = days,
                        sens_si_tot = sens_si_tot,
                        sens_si_nei = sens_si_nei,
                        sens_vax_tot = sens_vax_tot,
                        sens_vax_nei = sens_vax_nei,
                        sens_incid_tot = sens_incid_tot,
                        sens_incid_nei = sens_incid_nei,
                        sim_loc = sim_loc,
                        sim_loc_worse = sim_loc_worse,
                        sim_loc_better = sim_loc_better,
                        sim_loc_nei = sim_loc_nei,
                        list_calib_tot = list_calib_tot,
                        list_calib_nei = list_calib_nei,
                        scores_agg = scores_agg
)
# Save the list
saveRDS(output_analysis, file = "Output/all_analysis.RDS")
