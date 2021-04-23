# Function to generate population centroids 
source("R/function_location_centroids.R")
# Functions to generate the covariate time series and the serial interval
source("R/function_coverage.R")
source("R/function_distance_population.R")
source("R/function_serial_interval.R")
# Functions to simulate the outbreaks
source("R/function_generate_outbreak.R")
source("R/function_generate_all_outbreaks.R")
# Functions to analyse the outbreaks 
source("R/function_exp_gravity.R")
source("R/function_report_rate.R")
source("R/function_prepare_data_hhh4.R")
source("R/function_interpret_control_daily.R")
source("R/function_analysis_hhh4.R")
# Functions to generate the plots
source("R/function_figures.R")

# Import data table linking region nb and region id 
corres <- data.table(read.csv2(file = "Data/nuts_to_dep.csv", sep = ",", 
                               dec = ".", stringsAsFactors = F)[,c(2,1)])
colnames(corres) <- c("GeoCode", "regions")
setkey(corres, regions)
# Create covariate time series
cov_ts <- importation_coverage(corres = corres)
regs <- colnames(cov_ts)
all_dates <- as.Date(rownames(cov_ts))
# Create distance matrix using population centroids
dist_mat_cent <- importation_distance(regs = regs, neighbours = F)
# Create distance matrix using neighbours
dist_mat_nei <- importation_distance(regs = regs, neighbours = T)
# Create population and area time series
list_pop_area <- importation_pop_area(corres = corres, regs = regs, 
                                      all_dates = all_dates)
pop_ts <- list_pop_area$pop
area_ts <- list_pop_area$area

# Generate serial interval
w_dens <- serial_interval(prop_gen1 = 1, mean_w1 = 11.7, sd_w1 = 2.0)
