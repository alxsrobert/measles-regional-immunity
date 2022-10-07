importation_coverage <- function(corres, len_break = 1, rd_cov = FALSE, 
                                 quant_cov = 0.5,
                                 reg_exc = c("FRA10", "FRA20", "FRA30", "FRA40", 
                                             "FR831", "FR832", "FR")){
  # Import coverage csv file
  cov_france <- as.data.table(read.csv2(file = paste0(
    "Data/Vacc_coverage_departement.csv"),
    header = TRUE, stringsAsFactors = FALSE, sep = ",", dec = "."))
  # 2009 is missing
  colnames(cov_france) <- c("regions",
                            paste0("Year_", 2004:2008, "_1st"),
                            paste0("Year_", rep(2010:2017, each = 2), c("_1st", "_2nd")))
  # Extract values of 1st dose coverage
  table_coverage <- cov_france[, c(1, grep("_1st", colnames(cov_france))), with = FALSE]
  colnames(table_coverage) <- c("regions",
                                paste0("Year_", 2004:2008),
                                paste0("Year_", 2010:2017))

  if(class(table_coverage$regions) == "factor"){
    table_coverage[, regions := as.character(levels(regions)[regions])]
  }
  if(class(table_coverage$regions) != "character"){
    table_coverage[, regions := as.character(regions)]
  }
  # Convert all entries to numeric format
  table_coverage <- table_coverage[, lapply(.SD, as.numeric), by = regions]
  
  # Add column with region ID
  table_coverage[, Geocode := corres[table_coverage$regions, GeoCode]]
  # Rename column names
  years_col <- colnames(table_coverage)[grep("Year", colnames(table_coverage))]
  years_cov <- gsub(pattern = "Year_", years_col, replacement = "")
  colnames(table_coverage) <- c("regions", years_cov, "GeoCode")
  
  # Define starting and ending point of the dataset
  start_day <- as.Date(paste0(min(years_cov), "-01-01"))
  end_day <- as.Date(paste0(max(years_cov), "-12-31"))
  
  # Generate all dates
  break_days <- seq(start_day, end_day + 1, len_break)
  all_weeks <- unique(cut(as.Date(as.Date(start_day):as.Date(end_day),
                                  origin = "1970-01-01"),
                          breaks = break_days)) %>% as.character
  nweek <- length(all_weeks)
  
  # Data table daily coverage
  week_cov_dens <- data.table(week = rep(all_weeks, each = nrow(corres)),
                              regions = rep(table_coverage$GeoCode, nweek),
                              stringsAsFactors = FALSE)
  setkey(week_cov_dens, week)
  
  # Infer missing values
  dt_cov <- data.table(year = rep(as.numeric(years_cov), 
                                  each = nrow(table_coverage)), 
                       region = rep(table_coverage$GeoCode, length(years_cov)),
                       cov = unlist(c(table_coverage[, years_cov, with = FALSE])))
  dt_cov <- inference_coverage(corres = corres, dt_cov = dt_cov, 
                               rd_cov = rd_cov, quant_cov = quant_cov, 
                               reg_exc = reg_exc)
  
  ## Generate 3 year average
  # Create ID in dt_cov to match dt_cov to week_cov_dens
  dt_cov[, ID := paste0(region, "_", year)]
  setkey(dt_cov, ID)
  year_dt_week <- lubridate::year(as.Date(week_cov_dens$week))
  
  # For each entry, compute the local value of coverage 1, 2, and 3 years ago
  week_cov_dens[, cov_1 := dt_cov[paste0(regions, "_", year_dt_week - 1), cov]]
  week_cov_dens[, cov_2 := dt_cov[paste0(regions, "_", year_dt_week - 2), cov]]
  week_cov_dens[, cov_3 := dt_cov[paste0(regions, "_", year_dt_week - 3), cov]]
  
  # For each entry, compute the number of values of coverage in the past 3 years
  # that were reported or inferred, to take into account that no value was 
  # reported in 2009
  week_cov_dens[, n_cov := 3 - is.na(cov_1) - is.na(cov_2) - is.na(cov_3)]
  
  # Compute the average coverage
  week_cov_dens[is.na(cov_1), cov_1 := 0]
  week_cov_dens[is.na(cov_2), cov_2 := 0]
  week_cov_dens[is.na(cov_3), cov_3 := 0]
  week_cov_dens[, coverage := (cov_1 + cov_2 + cov_3)/n_cov]
  week_cov_dens <- week_cov_dens[!is.na(coverage),]
  # Remove the columns cov_1, cov_2, cov_3, and n_cov
  week_cov_dens[, c("cov_1", "cov_2", "cov_3", "n_cov") := NULL]

  regs <- unique(week_cov_dens$regions)
  cov_ts <- sapply(regs[is.na(match(regs, reg_exc))], function(X)
    return(week_cov_dens[regions == X, coverage]))
  rownames(cov_ts) <- week_cov_dens$week %>% unique %>% as.character
  colnames(cov_ts) <- regs[is.na(match(regs, reg_exc))]
  return(cov_ts)
}

inference_coverage <- function(corres, dt_cov, rd_cov = TRUE, quant_cov = 0.5,
                               reg_exc = c("FRA10", "FRA20", "FRA30", "FRA40", 
                                           "FR831", "FR832", "FR")){
  # Compute local coverage as a proportion
  dt_cov[, cov_prop := (cov/100)]
  # Set "year0", the number of year since 2003 
  dt_cov[, year0 := year - min(year) + 1]
  
  # Compute orthogonal polynomials of degree 2
  year_max <- max(dt_cov$year0)
  poly_year <- poly(seq_len(year_max), degree = 2)[dt_cov$year0,]
  
  # Add the orthogonal poly to dt_cov
  dt_cov[, year0_2p := poly_year[,2]]
  dt_cov[, year0_1p := poly_year[,1]]
  
  # Fit beta mixed model to the local proportion of vaccine coverage
  lmer2 <- glmmTMB(data = dt_cov[!is.element(region, reg_exc), ], 
                   family = beta_family(), 
                   cov_prop ~ 1 + year0_1p + year0_2p + 
                     (1 + year0_1p + year0_2p | region), se = TRUE)

  # Predict the mean values of the model for the missing entries
  preds <- predict(lmer2, se.fit = TRUE,
                   newdata = dt_cov[!is.element(region, reg_exc) & is.na(cov),
                                    .(region, year0_1p, year0_2p)], 
                   type = "conditional")
  # if rd_cov is true, draw the missing values from a normal distribution
  # Otherwise, use the mean value
  if(rd_cov){
    test_val <- rnorm(n = length(preds$fit), mean = preds$fit, sd = preds$se.fit) 
  } else test_val <- qnorm(p = quant_cov, mean = preds$fit, sd = preds$se.fit)
  
  # Add inferred values to dt_cov
  dt_cov[!is.element(region, reg_exc) & is.na(cov), cov := test_val * 100]
  
  # Remove columns year0 and cov_prop
  dt_cov[, c("year0", "cov_prop") := NULL]
  
  
  return(dt_cov)
  
  
}