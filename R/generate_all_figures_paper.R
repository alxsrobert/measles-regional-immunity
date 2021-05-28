########################################################
#### Generate all figures Main + figures Supplement ####
########################################################

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
# Import the analysis ran in generate_analysis_paper.R
list_analysis <- readRDS("Output/all_analysis.RDS")

#### Extract elements from list_analysis ####
# Reported case counts
data <- list_analysis$data_rep
# Model fits (1 and 2)
model_tot <- list_analysis$model_tot
model_nei <- list_analysis$model_nei
# Model fits in the sensitivity analysis (day of the week effect)
model_tot_weekday <- list_analysis$model_tot_weekday
# Simulations 1 year ahead (Models 1 and 2)
sim_1y_tot <- list_analysis$sim_1y_tot
sim_1y_nei <- list_analysis$sim_1y_nei
# Time span of the 1-year ahead simulations
days <- list_analysis$days
# Sensitivity analysis: impact of the serial interval (Models 1 and 2)
sens_si_tot <- list_analysis$sens_si_tot
sens_si_nei <- list_analysis$sens_si_nei
# Sensitivity analysis: impact of the inferred vaccine uptake (Models 1 and 2)
sens_vax_tot <- list_analysis$sens_vax_tot
sens_vax_nei <- list_analysis$sens_vax_nei
# Simulations 1-year ahead with group importations basic scenario
sim_loc <- list_analysis$sim_loc
# Simulations 1-year ahead with group importations variations in vaccine cov
sim_loc_worse <- list_analysis$sim_loc_worse
sim_loc_better <- list_analysis$sim_loc_better
# Simulations 1-year ahead with group importations model 2
sim_loc_nei <- list_analysis$sim_loc_nei
# Calibration scores Model 1, Model 2 and aggregated model 
list_calib_tot <- list_analysis$list_calib_tot
list_calib_nei <- list_analysis$list_calib_nei
scores_agg <- list_analysis$scores_agg

#### Generate figures Main ####

map <- generate_map()

# Figure 1: Data description
figure_1(data = data, mean = 11.7, sd = 2.0)

# Figure 2: Parameter estimates
figure_2(tot = model_tot, nei = model_nei, 
         ar_labs = c("Unvax", "Incid 1", 
                     "Incid 2", "Pop", 
                     "Area", "Sinus", "Cosinus"),
         ne_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                     "Sinus", "Cosinus"), 
         end_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                      "Sinus", "Cosinus"),
         other_labs = c("Distance", "Population", "Overdispersion"))

# Figure 3: Geographical distribution of the predictors
figure_3(model1 = model_tot, model2 = model_nei, map = map)

# Figure 4: Fit and calibration (Model 1)
figure_4(model = model_tot, list_calib = list_calib_tot)

# Figure 5: 1-year-ahead simulations: impact of coverage and incidence (Model 1)
# Create breaks, and labels for the figures
breaks <- c(-1, 1, 10, 50, 75, 90, 100)
labs <- c("0-1", "1-10", "10-50", "50-90", "90-99", "99-100")
# Create the label of each simulation set
name_elements <- c("Coverage 2018",
                   "Coverage 2018 plus 3%",
                   "Coverage 2018 minus 3%",
                   "Coverage 2018 &\nNo case in the last 3 years")

thresh_cases <- c(1, 10, 50)
order_reg = c(4, 3, 1, 2)

p <- figure_5_6(list_sim = sim_1y_tot, model = model_tot, map = map, 
                which_dates = seq_along(days), breaks = breaks, labs = labs,
                name_elements = name_elements, thresh_cases = thresh_cases, 
                order_reg = order_reg)
ggsave("Output/Figure_5.png", plot = p, width = 9, height = 12, units = "in")

# Figure 6: 1-year-ahead simulations: impact of group importations (Model 1)
regions <- c("FR623", "FR624", "FR824", "FR101")
types <- c("Haute-Garonne", "Gers", "Bouches-du-Rhône", "Paris")
import_loc <- coordinates(regions, types)
p <- figure_5_6(list_sim = sim_loc, model = model_tot, map = map, 
                which_dates = seq_along(days), breaks = breaks, labs = labs,
                name_elements = types, thresh_cases = thresh_cases, 
                order_reg = order_reg, import_loc = import_loc)
ggsave("Output/Figure_6.png", plot = p, width = 9, height = 12, units = "in")

#### Generate figures Supplement ####

# Figure S1: Sensitivity of the param estimates to different serial intervals
# (Model 1)
figures_sensitivity(model = model_tot, mat_coef = sens_si_tot$coefs,
                    min_coef = sens_si_tot$min_coef, 
                    max_coef = sens_si_tot$max_coef,
                    ar_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Area", "Sinus", "Cosinus"),
                    ne_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Sinus", "Cosinus"), 
                    end_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                 "Sinus", "Cosinus"),
                    other_labs = c("Distance", "Population", "Overdispersion"))

# Figure S2: Sensitivity of the param estimates to different serial intervals
# (Model 1) 
figures_sensitivity(model = model_nei, mat_coef = sens_si_nei$coefs,
                    min_coef = sens_si_nei$min_coef, 
                    max_coef = sens_si_nei$max_coef,
                    ar_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Area", "Sinus", "Cosinus"),
                    ne_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Sinus", "Cosinus"), 
                    end_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                 "Sinus", "Cosinus"),
                    other_labs = c("Population", "Overdispersion"))

# Figures S3, S4, S5, and S6: Figures describing the coverage beta mixed model
figures_coverage(corres = corres)

# Figure S7: Sensitivity of the param estimates to different vaccine uptakes
# (Model 1)
figures_sensitivity(model = model_tot, mat_coef = sens_vax_tot$coefs,
                    min_coef = sens_vax_tot$min_coef, 
                    max_coef = sens_vax_tot$max_coef,
                    ar_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Area", "Sinus", "Cosinus"),
                    ne_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Sinus", "Cosinus"), 
                    end_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                 "Sinus", "Cosinus"),
                    other_labs = c("Distance", "Population", "Overdispersion"))

# Figure S8: Sensitivity of the param estimates to different vaccine uptakes
# (Model 2)
figures_sensitivity(model = model_nei, mat_coef = sens_vax_nei$coefs,
                    min_coef = sens_vax_nei$min_coef, 
                    max_coef = sens_vax_nei$max_coef,
                    ar_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Area", "Sinus", "Cosinus"),
                    ne_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                "Sinus", "Cosinus"), 
                    end_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", 
                                 "Sinus", "Cosinus"),
                    other_labs = c("Population", "Overdispersion"))

# Figure S9: Seasonality in both models
figure_s9(model_tot, model_nei)

# Figure S10: Fit and calibration (Model 2)
figure_4(model = model_nei, list_calib = list_calib_nei)

# Figure S11: Trajectories of the calibration simulations
figure_s11(list_calib_tot = list_calib_tot, list_calib_nei = list_calib_nei,
           model = model_tot)

# Figure S12: Proportion per component
figure_s12(model_tot, model_nei)

# Figure S13: 1-year-ahead simulations: impact of coverage and incidence (Model 2)
name_elements <- c("Coverage 2018",
                   "Coverage 2018 plus 3%",
                   "Coverage 2018 minus 3%",
                   "Coverage 2018 &\nNo case in the last 3 years")
p <- figure_5_6(list_sim = sim_1y_nei, model = model_nei, map = map, 
                which_dates = seq_along(days), breaks = breaks, labs = labs,
                name_elements = name_elements, thresh_cases = thresh_cases, 
                order_reg = order_reg)
ggsave("Output/Figure_S13.png", plot = p, width = 9, height = 12, units = "in")

# Figure S14: 1-year-ahead simulations: impact of group importations (Model 2)
p <- figure_5_6(list_sim = sim_loc_nei, model = model_nei, map = map, 
                which_dates = seq_along(days), breaks = breaks, labs = labs,
                name_elements = types, thresh_cases = thresh_cases, 
                order_reg = order_reg, import_loc = import_loc)
ggsave("Output/Figure_S14.png", plot = p, width = 9, height = 12, units = "in")

# Figure S15: Last values of the covariates
figure_s15(map = map, 
           pop_mat = model_tot$control$data$pop, 
           unvax_mat = model_tot$control$data$unvax, 
           cat_incidence1 = model_tot$control$data$cat1,
           cat_incidence2 = model_tot$control$data$cat2)

# Figure S16: 1-year-ahead simulations: impact of group importations (Lower coverage)
p <- figure_5_6(list_sim = sim_loc_worse, model = model_tot, map = map, 
                which_dates = seq_along(days), breaks = breaks, labs = labs,
                name_elements = types, thresh_cases = thresh_cases, 
                order_reg = order_reg, import_loc = import_loc)
ggsave("Output/Figure_S16.png", plot = p, width = 9, height = 12, units = "in")

# Figure S17: 1-year-ahead simulations: impact of group importations (Higher coverage)
p <- figure_5_6(list_sim = sim_loc_better, model = model_tot, map = map, 
                which_dates = seq_along(days), breaks = breaks, labs = labs,
                name_elements = types, thresh_cases = thresh_cases, 
                order_reg = order_reg, import_loc = import_loc)
ggsave("Output/Figure_S17.png", plot = p, width = 9, height = 12, units = "in")

# Figure S18: Comparison aggregated and daily scores
figure_s18(scores_agg = scores_agg, scores_day = list_calib_tot[[3]])

# Figure S19: Number of cases per day of the weeks
figure_s19(data)

# Figure S20: Sensitivity of the param estimates to including day of the week effect
figure_2(model_tot_weekday, model_tot,
         ar_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", "Area", "Weekday", 
                     "Sinus", "Cosinus"),
         ne_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", "Weekday", "Sinus", 
                     "Cosinus"), 
         end_labs = c("Unvax", "Incid 1", "Incid 2", "Pop", "Weekday","Sinus", 
                      "Cosinus"),
         other_labs = c("Distance", "Population", "Overdispersion"),
         legend_text = c("Model 1 + Weekday effect", "Model 1"))
