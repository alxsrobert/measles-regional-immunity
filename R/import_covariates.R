source("R/function_location_centroids.R")

source("R/function_coverage.R")
source("R/function_distance_population.R")
source("R/function_serial_interval.R")

source("R/function_generate_outbreak.R")
source("R/function_generate_all_outbreaks.R")


corres <- data.table(read.csv2(file = "Data/nuts_to_dep.csv", sep = ",", 
                               dec = ".", stringsAsFactors = F)[,c(2,1)])

colnames(corres) <- c("GeoCode", "regions")
setkey(corres, regions)

cov_ts <- importation_coverage(corres = corres)
regs <- colnames(cov_ts)
all_dates <- as.Date(rownames(cov_ts))

dist_mat_cent <- importation_distance(regs = regs, neighbours = F)
list_pop_area <- importation_pop_area(corres = corres, regs = regs, 
                                      all_dates = all_dates)
pop_ts <- list_pop_area$pop
area_ts <- list_pop_area$area
