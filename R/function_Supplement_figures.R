figures_coverage <- function(corres, len_break = 1, 
                             reg_exc = c("FRA10", "FRA20", "FRA30", "FRA40", 
                                         "FR831", "FR832", "FR")){
  # Import coverage csv file
  cov_france <- as.data.table(read.csv2(file = paste0(
    "Data/Vacc_coverage_departement.csv"),
    header = T, stringsAsFactors = F, sep = ",", dec = "."))
  # 2009 is missing
  colnames(cov_france) <- c("regions",
                            paste0("Year_", 2004:2008, "_1st"),
                            paste0("Year_", rep(2010:2017, each = 2), c("_1st", "_2nd")))
  # Extract values of 1st dose coverage
  table_coverage <- cov_france[, c(1, grep("_1st", colnames(cov_france))), with = F]
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
                              stringsAsFactors = F)
  setkey(week_cov_dens, week)
  
  # Create a data table containing the year, region, and vaccine uptake
  dt_cov <- data.table(year = rep(as.numeric(years_cov), 
                                  each = nrow(table_coverage)), 
                       region = rep(table_coverage$GeoCode, length(years_cov)),
                       cov = unlist(c(table_coverage[, years_cov, with = F])))
  ## Generate S3: Show the number of missing entry per year / per region
  # Compute the number of missing entries per year
  missing_data <- dt_cov[!is.element(region, reg_exc), 
                         lapply(.SD, function(X) sum(is.na(X))), 
                         .SDcols = "cov", by = year]$cov
  names(missing_data) <- dt_cov[!is.element(region, reg_exc), 
                                lapply(.SD, function(X) sum(is.na(X))),
                                .SDcols = "cov", by = year]$year
  # Number of measurements per year
  year_data <- dt_cov[!is.element(region, reg_exc), year] %>% table
  # No measurement in 2009: Add 2009 to "missing_data", with 94 missing entries
  missing_data["2009"] <- unique(year_data)
  missing_data <- missing_data[order(names(missing_data))]
  # Compute the number of missing entries per region
  missing_region <- dt_cov[!is.element(region, reg_exc), 
                           lapply(.SD, function(X) sum(is.na(X))),
                           .SDcols = "cov", by = region]$cov
  
  par(mfrow = c(1,2), mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 1.5)
  ## Panel A: Number of missing entry per year
  plot(as.numeric(year_data)~as.numeric(names(year_data)), type = "h", 
       lwd = 10, lend = "butt", ylim = c(0, 110), xlab = "Year", 
       ylab = "Number of observations")
  lines(missing_data~as.numeric(names(missing_data)), type = "h", 
        lwd = 10, lend = "butt", col = "red")
  lines(missing_data[names(missing_data)<2006 | names(missing_data) == 2009]~
          as.numeric(names(missing_data[names(missing_data)<2006 | 
                                          names(missing_data) == 2009])), type = "h", 
        lwd = 10, lend = "butt", col = "grey")
  mtext("A", side = 3, line = -2, adj = 0.1, cex = 2)
  legend(lwd = 10, "topright", ncol = 1, legend = c("Reported", "Inferred", 
                                                    "Unreported, not inferred"),
         col = c("black", "red", "grey"), bty = "n")
  ## Panel B: Number of missing entry in each region
  plot(as.numeric(table(missing_region))~as.numeric(names(table(missing_region))),
       type = "p", pch = 16, xlab = "number of missing entries", 
       ylab = "Number of regions")
  mtext("B", side = 3, line = -2, adj = 0.1, cex = 2)
  
  ## Fit the mixed model to the coverage data
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
                     (1 + year0_1p + year0_2p | region), se = T)
  ## Plot the typical trajectories for 3 regions
  # Randomly draw the regions to plot
  set.seed(1)
  regions_for_plot <- sample(unique(dt_cov[!is.element(region, reg_exc),
                                           region]), 3)
  # Compute the trajectory across all regions
  pred_tot <- predict(lmer2, se.fit = T, re.form = NA,
                      newdata = dt_cov[!is.element(region, reg_exc) & 
                                         !duplicated(year),
                                       .(region, year0_1p, year0_2p)])
  # Compute mean prediction and 95% CIs
  mean_pred <- exp(pred_tot[[1]])/(1 + exp(pred_tot[[1]]))
  min_pred <- exp(qnorm(.025, mean = pred_tot[[1]], sd = pred_tot[[2]])) /
    (1 + exp(qnorm(.025, mean = pred_tot[[1]], sd = pred_tot[[2]])))
  max_pred <- exp(qnorm(.975, mean = pred_tot[[1]], sd = pred_tot[[2]])) /
    (1 + exp(qnorm(.975, mean = pred_tot[[1]], sd = pred_tot[[2]])))
  # Compute the trajectory for the three regions selected
  pred_reg <- predict(lmer2, se.fit = T,
                      newdata = dt_cov[is.element(region, regions_for_plot),
                                       .(region, year0_1p, year0_2p)])
  # Compute mean prediction and 95% CIs
  mean_pred_reg <- exp(pred_reg[[1]])/(1 + exp(pred_reg[[1]]))
  min_pred_reg <- exp(qnorm(.025, mean = pred_reg[[1]], sd = pred_reg[[2]])) /
    (1 + exp(qnorm(.025, mean = pred_reg[[1]], sd = pred_reg[[2]])))
  max_pred_reg <- exp(qnorm(.975, mean = pred_reg[[1]], sd = pred_reg[[2]])) /
    (1 + exp(qnorm(.975, mean = pred_reg[[1]], sd = pred_reg[[2]])))
  
  # Compute the years where at least 1 measure was reported 
  unique_year <- dt_cov[!is.element(region, reg_exc) & !duplicated(year), year]
  # Compute the years with reported data in the regions selected
  year_region <- dt_cov[is.element(region, regions_for_plot), year]
  
  
  # Generate the colour scheme (1 colour per region selected)
  cols <- brewer.pal(n = length(regions_for_plot), name = "Accent")
  names(cols) <- regions_for_plot
  
  reg_region <- dt_cov[is.element(region, regions_for_plot), region]
  
  par(mfrow = c(2,1), mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 1.5)
  # Panel A: Plot the values of the coeffients for each region
  pos <- order(coef(lmer2)[[1]][[1]][,1])
  cols_coef <- brewer.pal(n = 3, name = "Paired")
  plot(coef(lmer2)[[1]][[1]][pos,1], pch = 16, ylim = c(-.5, 4), 
       xlab = "Regions", ylab = "Coefficients", 
       col = cols_coef[1])
  points(coef(lmer2)[[1]][[1]][pos,2], pch = 16, col = cols_coef[2])
  points(coef(lmer2)[[1]][[1]][pos,3], pch = 16, col = cols_coef[3])
  legend("topleft", pch = 16, col = cols_coef, bty = "n", cex = 1.5,
         legend = c("Intercept", "Year orthogonal poly 1",
                    "Year orthogonal poly 2"))
  # Compute the 4 trajectories (Mean trajectory + 3 regions previously selected)
  plot(mean_pred~unique_year, pch = 18, ylim = c(.7, 1), ylab = "Coverage", 
       xlab = "Time")
  polygon(x = c(unique_year, rev(unique_year)), y = c(min_pred, rev(max_pred)), 
          col = transp("black", .3), border = NA)
  points(mean_pred_reg~year_region, pch = 18, col = cols[reg_region])
  for(i in regions_for_plot){
    pos <- which(reg_region == i)
    polygon(x = c(unique_year, rev(unique_year)), y = c(min_pred_reg[pos],
                                                        rev(max_pred_reg[pos])), 
            col = transp(cols[i], .3), border = NA)
  }
  legend("topleft", pch = 16, col = c("black", cols), 
         legend = c("Average", "Region 1", "Region 2", "Region 3"),
         bty = "n", cex = 1.5, ncol = length(regions_for_plot) + 1)
  
  points(mean_pred~unique_year, pch = 18, ylim = c(.80, 1))
  polygon(x = c(unique_year, rev(unique_year)), y = c(min_pred, rev(max_pred)), 
          col = transp("black", .3), border = NA)
  
  ## Generate diagnostic plots of the mixed-effect regression
  # Simulate values using the coefficients from the model
  lmer2_res <- simulateResiduals(lmer2, plot = F, quantreg = F, n_sim = 1000,
                                 seed = 1)
  par(mfrow = c(1,2), mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 1.5)
  # Panel A: Fitted vs Residuals
  plot(residuals(lmer2)~fitted(lmer2), pch = 16, col = transp("black", 0.3),
       ylim = c(-0.1, 0.1), ylab = "Residuals", xlab = "Fitted")
  abline(h = 0, col = "red", lty = 2)
  # Panel B: QQ plots
  plotQQunif(lmer2_res, testUniformity = F, testOutliers = F, testDispersion = F)
  
  # dd: Full dataset
  dd <- dt_cov[!is.na(cov) & !is.element(region, reg_exc),
               .(year0_1p, year0_2p, region, cov, year)]
  # aa: Fitted values and residuals for each data point 
  aa <- augment(lmer2, data=dd)
  # fitted0 = fitted excluding the random effects
  aa$.fitted0 <- predict(lmer2, newdata=transform(dd,region=NA),type="response")
  # resid0 = Residuals excluding the random effects
  aa$.resid0 <- dd$cov/100-aa$.fitted0
  # Plot the residuals for each year
  gg3 <- (ggplot(aa, aes(year,.resid0))
          + geom_line(aes(group=region),colour= transp("lightgrey", .3))
          + geom_point() + geom_smooth() + ylim(-.2, .2) 
          + theme_classic(base_size = 9) + labs(x = "Year", y = "Residuals")
  )
  plot(gg3)
}  

figure_s9 <- function(tot, nei){
  ## Extract the two seasonality coefficients for each component:
  #-in tot
  tot_cos_ar <- coef(tot)["ar.cos(2 * pi * t/365)"]
  tot_cos_ne <- coef(tot)["ne.cos(2 * pi * t/365)"]
  tot_cos_en <- coef(tot)["end.cos(2 * pi * t/365)"]
  tot_sin_ar <- coef(tot)["ar.sin(2 * pi * t/365)"]
  tot_sin_ne <- coef(tot)["ne.sin(2 * pi * t/365)"]
  tot_sin_en <- coef(tot)["end.sin(2 * pi * t/365)"]
  #- in nei
  nei_cos_ar <- coef(nei)["ar.cos(2 * pi * t/365)"]
  nei_cos_ne <- coef(nei)["ne.cos(2 * pi * t/365)"]
  nei_cos_en <- coef(nei)["end.cos(2 * pi * t/365)"]
  nei_sin_ar <- coef(nei)["ar.sin(2 * pi * t/365)"]
  nei_sin_ne <- coef(nei)["ne.sin(2 * pi * t/365)"]
  nei_sin_en <- coef(nei)["end.sin(2 * pi * t/365)"]
  t <- (as.Date("2019-01-01") + seq_len(366) - 1)
  t_num <- (t - tot$control$data$unvax %>% rownames %>% as.Date %>% min) %>% 
    as.numeric
  
  par(mfrow = c(2,1), mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 1.5)
  plot(100*(exp(tot_cos_ar * cos(2 * pi * t_num/365) + 
                  tot_sin_ar * sin(2 * pi * t_num/365)) - 1) ~ 
         t, type = "l", col = "orange", lwd = 2, 
       xlab = "Time", ylab = "Percent", ylim = c(-55, 55))
  lines(100*(exp(tot_cos_ne * cos(2 * pi * t_num/365) + 
                   tot_sin_ne * sin(2 * pi * t_num/365))-1) ~ 
          t, col = "blue", lwd = 2)
  lines(100*(exp(tot_cos_en * cos(2 * pi * t_num/365) + 
                   tot_sin_en * sin(2 * pi * t_num/365))-1) ~ 
          t, col = "grey", lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
  mtext("A", side = 3, line = -1, adj = 0.1, cex = 2)
  legend("bottomleft", col = c("orange", "blue", "grey"), lwd = 2,
         legend = c("Autoregressive", "Neighbour", "Endemic"), bty = "n")
  plot(100*(exp(nei_cos_ar * cos(2 * pi * t_num/365) + 
                  nei_sin_ar * sin(2 * pi * t_num/365)) - 1) ~ 
         t, type = "l", col = "orange", lwd = 2, 
       xlab = "Time", ylab = "Percent", ylim = c(-55, 55))
  lines(100*(exp(nei_cos_ne * cos(2 * pi * t_num/365) + 
                   nei_sin_ne * sin(2 * pi * t_num/365))-1) ~ 
          t, col = "blue", lwd = 2)
  lines(100*(exp(nei_cos_en * cos(2 * pi * t_num/365) + 
                   nei_sin_en * sin(2 * pi * t_num/365))-1) ~ 
          t, col = "grey", lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
  mtext("B", side = 3, line = -1, adj = 0.1, cex = 2)
  
}

figure_s11 <- function(tot, nei){
  tot_hhh <- meanHHH(coef(tot), terms(tot))
  tot_nb <- c(Endemic = sum(tot_hhh$endemic), 
              Autoregressive = sum(tot_hhh$epi.own),
              Neighbour = sum(tot_hhh$epi.neighbours))
  tot_prop <- tot_nb/sum(tot_nb) * 100
  nei_hhh <- meanHHH(coef(nei), terms(nei))
  nei_nb <- c(Endemic = sum(nei_hhh$endemic), 
              Autoregressive = sum(nei_hhh$epi.own),
              Neighbour = sum(nei_hhh$epi.neighbours))
  nei_prop <- nei_nb/sum(nei_nb) * 100
  
  
  par(mfrow = c(1,1), mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 2)
  barplot(rbind(tot_prop, nei_prop), beside = T, col = c("blue", "purple"),
          main = "Proportion of cases per component", ylab = "Percentage",
          ylim = c(0, 100))
  legend("topright", legend = c("Model 1: Spatial spread between all regions", 
                                "Model 2: Spatial spread between neighbours only"),
         lwd = 2, col = c("blue", "purple"), bty = "n", cex  = 1.1)
  
}

figure_s14 <- function(pop_mat, unvax_mat, cat_incidence1, cat_incidence2, map){
  map$pop <- pop_mat[nrow(pop_mat), ][map$nuts3]
  map$pop_cat <- cut(map$pop, breaks = c(0, 300000, 500000, 750000, 1000000, 
                                         max(map$pop, na.rm = T)), 
                     labels = c("70k-300k", "300k-500k", "500k-750k", "750k-1M", ">1M"))
  p_cases1 <- ggplot(map)  +  geom_sf(aes(fill = pop_cat)) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
                      name = "Number of inhabitants") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), legend.position = "bottom", 
          axis.line = element_blank()) + guides(fill = guide_legend(ncol = 2))
  map$cov <- 100 * (1 - exp(unvax_mat[nrow(unvax_mat), ]))[map$nuts3]
  map$cov_cat <- cut(map$cov, breaks = c(0, 80, 85, 90, 95, 
                                         max(map$cov, na.rm = T)), 
                     labels = c("75%-80%", "80%-85%", "85%-90%", "90%-95%", ">95%"))
  p_cases2 <- ggplot(map)  +  geom_sf(aes(fill = cov_cat)) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
                      name = "Vaccine coverage") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), legend.position = "bottom", 
          axis.line = element_blank()) + guides(fill = guide_legend(ncol = 2))
  
  map$incidence <- cat_incidence1[nrow(pop_mat), ][map$nuts3] + 
    2*cat_incidence2[nrow(pop_mat), ][map$nuts3]
  map$incidence_cat <- cut(map$incidence, breaks = c(-1, 0, 1, 
                                                     max(map$incidence, na.rm = T)), 
                           labels = c("No recent outbreak", 
                                      "Moderate recent transmission",
                                      "High recent transmission"))
  p_cases3 <- ggplot(map)  +  geom_sf(aes(fill = incidence_cat)) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#ffffbf", "#fdae61", "#d7191c"),
                      name = "Recent incidence") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), legend.position = "bottom", 
          axis.line = element_blank()) + guides(fill = guide_legend(ncol = 1))
  test_arrange <- arrangeGrob(p_cases1, p_cases2, p_cases3, ncol = 3)
  plot(test_arrange)
}