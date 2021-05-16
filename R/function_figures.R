# Plots describing the simulations
plot_simulations <- function(list_all_sim, list_all_sim_reported, map, 
                             breaks_vec = c(-1, 5, 25, 50, 75, 95, 100),
                             labs_vec = c("0-5", "6-25", "26-50", "51-75", 
                                          "76-95", "96-100"),
                             thresh = c(1, 10, 50), which_plots = c(1,2)){
  # Compute the overall number of cases generated
  n_per_sim <- list_all_sim[[1]] %>% lapply(sum) %>% unlist %>% sort
  # Compute the number of cases reported
  n_per_sim_rep <- list_all_sim_reported[[1]] %>% lapply(sum) %>% unlist %>% sort
  # Vector showing the year of each date entry in the simulation
  year_sim <- rownames(list_all_sim_reported[[1]][[1]]) %>% as.Date %>% 
    lubridate::year()
  # Compute the number of cases generated per year
  n_per_year <- list_all_sim[[1]] %>% lapply(function(X){
    n_per_year <- aggregate(x = rowSums(X), by = list(year_sim), sum)$x
    names(n_per_year) <- unique(year_sim)
    return(n_per_year)
  }) %>% do.call(rbind,.)
  # Compute the number of cases reported per year
  n_per_year_rep <- list_all_sim_reported[[1]] %>% lapply(function(X){
    n_per_year <- aggregate(x = rowSums(X), by = list(year_sim), sum)$x
    names(n_per_year) <- unique(year_sim)
    return(n_per_year)
  }) %>% do.call(rbind,.)
  
  
  rep_per_region_1 <- list_all_sim_reported[[1]] %>% lapply(function(X){
    ab_1 <- (aggregate(x = X, by = list(year_sim), sum)[,-1] >= thresh[1]) %>% 
      apply(2, any)
    return(ab_1)
  }) %>% do.call(rbind,.) %>% colSums()
  rep_per_region_2 <- list_all_sim_reported[[1]] %>% lapply(function(X){
    ab_2 <- (aggregate(x = X, by = list(year_sim), sum)[,-1] >= thresh[2]) %>% 
      apply(2, any)
    return(ab_2)
  }) %>% do.call(rbind,.) %>% colSums()
  rep_per_region_3 <- list_all_sim_reported[[1]] %>% lapply(function(X){
    ab_3 <- (aggregate(x = X, by = list(year_sim), sum)[,-1] >= thresh[3]) %>% 
      apply(2, any)
    return(ab_3)
  }) %>% do.call(rbind,.) %>% colSums()
  
  map1 <- map2 <- map3 <- map
  map1$ab <- rep_per_region_1[as.character(map$nuts)]
  map2$ab <- rep_per_region_2[as.character(map$nuts)]
  map3$ab <- rep_per_region_3[as.character(map$nuts)]
  if(thresh[1] == 1) map1$type <- paste0("At least 1 case") else
    map1$type <- paste0("At least ", thresh[1], " cases")
  map2$type <- paste0("At least ", thresh[2], " cases")
  map3$type <- paste0("At least ", thresh[3], " cases")
  
  # Cut in categories
  map1$ab_cat <- cut(map1$ab, breaks = breaks_vec, labels = labs_vec)
  map2$ab_cat <- cut(map2$ab, breaks = breaks_vec, labels = labs_vec)
  map3$ab_cat <- cut(map3$ab, breaks = breaks_vec, labels = labs_vec)
  maptot <- rbind(map1, map2, map3)
  
  par(mfrow = c(2,2), mar = c(5,5,1,1), las = 1,
      bty ="l", cex.axis = 1.1, cex.main = 1.5, cex.lab = 1.5)
  # Generate number of cases in each simulation set
  plot(n_per_sim, pch = 16, log = "y", ylim = c(3000, 300000), 
       xlab = "Simulation", ylab = "", col = transp("black", .5))
  title(ylab = "Number of cases", line = 4)
  points(n_per_sim_rep, pch = 16, , col = transp("red", .5))
  legend("left", legend = c("Cases generated", "Cases reported"),
         border = NA, fill = c(transp("black", .7), transp("red", .7)),
         bty = "n", cex = 1.3)
  mtext("A", side = 3, line = -2, adj = 0.1, cex = 2)
  # Generate boxplots of the number of cases per year
  boxplot(n_per_year_rep, log = "y", xlab = "Year", pch = 16, 
          col = transp("red", .5))
  title(ylab = "Number of cases", line = 4)
  mtext("B", side = 3, line = -2, adj = 0.1, cex = 2)
  time <- as.Date(rownames(list_all_sim[[1]][[1]]))
  # Generate time series for two of the simulation sets
  plot(rowSums(list_all_sim[[1]][[which_plots[1]]]) ~ time, type = "l", 
       ylab = "Number of cases", xlab = "Day", col = transp("black", .5))
  lines(rowSums(list_all_sim_reported[[1]][[which_plots[1]]]) ~ time, 
        col = transp("red", .5))
  mtext("C", side = 3, line = -2, adj = 0.1, cex = 2)
  plot(rowSums(list_all_sim[[1]][[which_plots[2]]]) ~ time, type = "l",
       ylab = "Number of cases", xlab = "Day", col = transp("black", .5))
  lines(rowSums(list_all_sim_reported[[1]][[which_plots[2]]]) ~ time, 
        col = transp("red", .5))
  mtext("D", side = 3, line = -2, adj = 0.1, cex = 2)
  # Generate number of cases per region in the simulations
  p_cases <- ggplot(maptot) +  geom_sf(aes(fill = ab_cat)) + 
    facet_grid(.~type) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#2c7bb6", "#abd9e9", "#e0f3f8", 
                                 "#ffffbf", "#fdae61", "#d7191c"),
                      name = "Percentage of simulations \n with one year above threshold") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), strip.text = element_text(size = 12),  
          axis.ticks = element_blank(), axis.line = element_blank(), 
          legend.position = "bottom") + guides(fill = guide_legend(byrow = T))
  return(p_cases)
}
# Function to add transparency to colours 
transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) 
    rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}
# Function to compute the difference mean estimates - parameters
distance_params <- function(list_all_sim, hhh4_runs){
  # Extract the parameters used to generate the simulations
  params_mat <- list_all_sim$params_sim
  # Initialise difference matrix 
  dist <- t(params_mat) * 0
  
  for (i in seq_along(hhh4_runs)){
    # Extract coefficient estimated in run i
    coef_i <- coef(hhh4_runs[[i]])
    # Make sure overdisp is on the same scale as params_mat
    coef_i["overdisp"] <- 1/coef_i["overdisp"]
    # Compute the difference
    dist[,i] <- as.numeric(coef_i - params_mat[i,])
  }
  return(dist)
}
# Function to compute whether the input parameters fall into the estimated CIs
in_CI <- function(list_all_sim, hhh4_runs, CI = 0.95){
  # Extract the parameters used to generate the simulations
  params_mat <- list_all_sim$params_sim
  
  # Initialise CI matrix 
  in_CI_mat <- t(params_mat) * 0
  
  
  for (i in seq_along(hhh4_runs)){
    # Extract CI% confidence intervals at run i
    ci_low <- confint(hhh4_runs[[i]], level = CI)[,1]
    ci_up <- confint(hhh4_runs[[i]], level = CI)[,2]
    # Make sure overdisp is on the same scale as params_mat
    up_overdisp <- 1/ci_low["overdisp"]
    low_overdisp <- 1/ci_up["overdisp"]
    ci_low["overdisp"] <- low_overdisp
    ci_up["overdisp"] <- up_overdisp
    # Binary value: 1 if the reference parameter is in the CI, 0 otherwise
    in_CI_mat[,i] <- params_mat[i,] > ci_low & params_mat[i,] < ci_up
  } 
  return(in_CI_mat)
}
# Function to compute the pit histogram of the parameters
pit_param <- function(list_all_sim, hhh4_runs){
  # Extract the parameters used to generate the simulations
  params_mat <- list_all_sim$params_sim
  
  # Initialise pit matrix 
  pit_mat <- t(params_mat) * 0
  
  for (i in seq_along(hhh4_runs)){
    # Extract mean coefficients estimated in run i
    coefs_mod <- coef(hhh4_runs[[i]], reparamPsi = F)
    # Extract standard deviation of the coefficients estimated in run i
    sd_mod <- hhh4_runs[[i]]$se
    # Extract the parameter set at run i
    val_params <- list_all_sim$params_sim[i,]
    val_params["overdisp"] <- -log(1/val_params["overdisp"])
    # Compute the corresponding probability
    pit_mat[,i] <- pnorm(q = as.numeric(val_params), 
                         mean = coefs_mod, sd = sd_mod)
  } 
  return(pit_mat)
}
# Overall function to generate all figures
plot_analysis <- function(list_all_sim, hhh4_day, hhh4_agg, CI, which_plot, 
                          exclude_overdisp = T, mains = NULL){
  # If which_plot is not numeric, match the elements to the parameter names
  if(!is.numeric(which_plot)){
    which_plot <- match(which_plot, colnames(list_all_sim$params_sim))
  }
  # Generate CI matrices
  CI_day <- in_CI(list_all_sim = list_all_sim, hhh4_runs = hhh4_day, CI = CI)
  CI_agg <- in_CI(list_all_sim = list_all_sim, hhh4_runs = hhh4_agg, CI = CI)
  # Remove overdisp parameter
  if(exclude_overdisp){
    CI_day <- CI_day[-25,]
    CI_agg <- CI_agg[-25,]
  }
  # Remove intercepts
  CI_day <- CI_day[-grep("intercept", rownames(CI_day)),]
  CI_agg <- CI_agg[-grep("intercept", rownames(CI_agg)),]
  # Generate difference between mean of the estimates and the parameter set
  dist_day <- distance_params(list_all_sim = list_all_sim, hhh4_runs = hhh4_day)
  dist_agg <- distance_params(list_all_sim = list_all_sim, hhh4_runs = hhh4_agg)
  # Generate pit values
  pit_i_day <- pit_param(list_all_sim = list_all_sim, hhh4_runs = hhh4_day)
  pit_i_agg <- pit_param(list_all_sim = list_all_sim, hhh4_runs = hhh4_agg)
  # Histograms of the proportion of parameters within the confidence interval
  par(mfrow = c(2,1), mar = c(5,5,1,1), las = 1,
      bty ="l", cex.axis = 1.1, cex.main = 1.5, cex.lab = 1.5)
  hist(rowSums(CI_day)/ncol(CI_day) * 100, main = "", breaks = seq(0, 100, 2.5), 
       xlab = "", ylab = "")
  title(ylab = "Number of parameters")
  mtext("A: Daily models", side = 3, line = -2, adj = 0.1, cex = 2)
  abline(v = 95, col = "red", lwd = 2, lty = 2)
  hist(rowSums(CI_agg)/ncol(CI_agg) * 100, main = "", breaks = seq(0, 100, 2.5), 
       xlab = "", ylab = "")
  abline(v = 95, col = "red", lwd = 2, lty = 2)
  mtext("B: Aggregated models", side = 3, line = -2, adj = 0.1, cex = 2)
  title(ylab = "Number of parameters")
  title(xlab = 
          "Proportion of simulations with the input parameter included in the 95% confidence interval of the model")
  ## For each parameter listed in which_plot, generate 3 figures:
  # 1- The density plot of the difference between estimates and parameter 
  # 2- Density plots of the mean of the parameter (with data)
  # 3- pit histogram of the figure
  n_row <- min(c(3, length(which_plot)))
  par(mfrow = c(n_row, 3), mar = c(5, 5, 1, 1), las = 1,
      bty ="l", cex.axis = 1.1, cex.main = 1.5, cex.lab = 1.5)
  for(i in seq_along(which_plot)){
    if(i %% (n_row) == 1){
      par(mfrow = c(n_row, 3), mar = c(5, 5, 1, 1), las = 1, 
          bty ="l", cex.axis = 1.1, cex.main = 1.5, cex.lab = 1.5)
    }
    param_i <- which_plot[i]
    if(!is.null(mains)) main_i <- mains[i] else 
      main_i <- colnames(list_all_sim$params_sim)[param_i]
    # Generate density objects
    dens_dis_day <- density(dist_day[param_i,], na.rm = T)
    dens_dis_agg <- density(dist_agg[param_i,], na.rm = T)
    # Initialise empty plot
    plot(dens_dis_agg, col = "white", xlab = "", ylab = "", main = main_i, 
         ylim = c(0, max(c(dens_dis_agg$y, dens_dis_day$y))),
         xlim = c(max(-3, min(c(dens_dis_agg$x, dens_dis_day$x))),
                  max(c(dens_dis_agg$x, dens_dis_day$x))))
    if(i %% n_row == 0 || i == length(which_plot)) 
      title(xlab = "Distance to data")
    if(i %% n_row == 1) {
      legend("left", legend = c("Daily model", "Aggregated model", "Simulations"),
             border = NA, fill = c(transp("#fc8d59", .3), transp("#91bfdb", .3),
                                   transp("black", .3)),
             bty = "n", cex = 1.3)
      mtext("A", side = 3, line = -2, adj = 0.1, cex = 2)
    }
    if(i %% n_row == 2) mtext("D", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i %% n_row == 0) mtext("G", side = 3, line = -2, adj = 0.1, cex = 2)

    title(ylab = "Density")
    # Add density plots
    polygon(dens_dis_day, col=transp("#fc8d59", .3), border=NA)
    polygon(dens_dis_agg, col=transp("#91bfdb", .3), border=NA)
    abline(v = 0, lty = 2, lwd = 2)
    # Extract mean values of the coefficients 
    mean_day <- unlist(lapply(hhh4_day, function(X) return(coef(X)[param_i])))
    mean_agg <- unlist(lapply(hhh4_agg, function(X) return(coef(X)[param_i])))
    # Compute density
    dens_day <- density(mean_day, na.rm = T)
    dens_agg <- density(mean_agg, na.rm = T)
    # Compute density of the data 
    dens_data <- density(list_all_sim$params_sim[, param_i], na.rm = T)
    # Initialise empty plot
    plot(dens_agg, col = "white", xlab = "", ylab = "", main = main_i, 
         ylim = c(0, max(c(dens_day$y, dens_data$y, dens_agg$y))),
         xlim = c(max(c(min(c(dens_agg$x, dens_day$x, dens_data$x)), -3)), 
                  max(c(dens_agg$x, dens_day$x, dens_data$x))))
    # Add density plots
    polygon(dens_data, col=transp("black", .3), border=NA)
    polygon(dens_day, col=transp("#fc8d59", .3), border=NA)
    polygon(dens_agg, col=transp("#91bfdb", .3), border=NA)
    abline(v = 0, lty = 2, lwd = 2)
    if(i %% n_row == 0 || i == length(which_plot)) 
      title(xlab = "Parameter value")
    if(i %% n_row == 1) mtext("B", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i %% n_row == 2) mtext("E", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i %% n_row == 0) mtext("H", side = 3, line = -2, adj = 0.1, cex = 2)
    title(ylab = "Density")
    # Generate histograms
    hist_day <- hist(pit_i_day[param_i,], breaks = seq(0,1, .1), plot = F)
    hist_agg <- hist(pit_i_agg[param_i,], breaks = seq(0,1, .1), plot = F)
    hist_day$counts <- hist_day$density
    hist_agg$counts <- hist_agg$density
    plot(hist_day, col = transp("#fc8d59", .3), ylab = "", border = NA,
         main = main_i, xlab = "", 
         ylim = c(0, max(c(hist_day$counts, hist_agg$counts))))
    plot(hist_agg, col = transp("#91bfdb", .3), add = TRUE, border = NA)
    abline(h = 1, lwd = 2, lty = 2)
    if(i %% n_row == 0 || i == length(which_plot)) 
      title(xlab = "Probability Integral Transform")
    if(i %% n_row == 1) mtext("C", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i %% n_row == 2) mtext("F", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i %% n_row == 0) mtext("I", side = 3, line = -2, adj = 0.1, cex = 2)
    title(ylab = "Relative frequency")
  }
}
# Proportion of cases per component
plot_prop_comp <- function(list_all_sim, hhh4_day, hhh4_agg){
  # Extract the proportion of cases from each component in the daily model
  prop_daily <- lapply(models_daily, function(X){
    mean_X <- meanHHH(coef(X), terms(X))
    nb_comp <- c(Autoregressive = sum(mean_X$epi.own),
                 Neighbourhood = sum(mean_X$epi.neighbours),
                 Endemic = sum(mean_X$endemic))
    return(nb_comp / sum(nb_comp))
  }) %>% do.call(rbind, .)

  # Extract the proportion of cases from each component in the aggregated model
  prop_aggre <- lapply(models_aggre, function(X){
    mean_X <- meanHHH(coef(X), terms(X))
    nb_comp <- c(Autoregressive = sum(mean_X$epi.own),
                 Neighbourhood = sum(mean_X$epi.neighbours),
                 Endemic = sum(mean_X$endemic))
    return(nb_comp / sum(nb_comp))
  }) %>% do.call(rbind, .)
  # Generate the median and 95% CI for each component, daily model
  med_day <- apply(prop_daily, 2, median) * 100
  min_day <- apply(prop_daily, 2, function(X) quantile(X, probs = c(0.025))) * 100
  max_day <- apply(prop_daily, 2, function(X) quantile(X, probs = c(0.975))) * 100
  # Generate the median and 95% CI for each component, aggregated model
  med_agg <- apply(prop_aggre, 2, median) * 100
  min_agg <- apply(prop_aggre, 2, function(X) quantile(X, probs = c(0.025))) * 100
  max_agg <- apply(prop_aggre, 2, function(X) quantile(X, probs = c(0.975))) * 100
  # Generate the median and 95% CI for each component, data
  med_data <- apply(list_all_sim$prop_comp, 2, median) * 100
  min_data <- apply(list_all_sim$prop_comp, 2, 
                    function(X) quantile(X, probs = c(0.025))) * 100
  max_data <- apply(list_all_sim$prop_comp, 2, 
                    function(X) quantile(X, probs = c(0.975))) * 100
  
  # Use barplots to represent the medians, and arrows for the 95% CI
  par(mfrow = c(1,1), mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 2)
  b <- barplot(rbind(med_data, med_day, med_agg), beside = T, 
               col = c("grey", "blue", "purple"),
               main = "Proportion of cases per component", ylab = "Percentage",
               ylim = c(0, 100), border = NA)
  arrows(x0 = b[1,], y0 = min_data, x1 = b[1,], y1 = max_data, 
         angle = 90, code = 3, length = 0.05, lwd = 2)
  arrows(x0 = b[2,], y0 = min_day, x1 = b[2,], y1 = max_day, 
         angle = 90, code = 3, length = 0.05, lwd = 2)
  arrows(x0 = b[3,], y0 = min_agg, x1 = b[3,], y1 = max_agg, 
         angle = 90, code = 3, length = 0.05, lwd = 2)
  legend("topright", legend = c("Simulations", "Daily model", "Aggregated model"),
         fill = c("grey", "blue", "purple"), bty = "n", cex  = 1.1, border = NA)
  
}
# Difference in calibration scores
plot_calib_scores <- function(scores_agg, scores_day){
  # Generate mean value of bias per simulation, aggregated and daily models
  mean_bias_agg <- apply(scores_agg$bias, 3, mean)
  mean_bias_day <- apply(scores_day$bias, 3, mean)
  # Generate mean value of sharpness per simulation, aggregated and daily models
  mean_shar_agg <- apply(scores_agg$shar, 3, mean)
  mean_shar_day <- apply(scores_day$shar, 3, mean)
  # Generate mean value of RPS per simulation, aggregated and daily models
  mean_rps_agg <- apply(scores_agg$rps, 3, mean)
  mean_rps_day <- apply(scores_day$rps, 3, mean)
  ## Compute and categorise the difference between aggregated and daily models
  # Bias
  cat_diff_bias <- cut(mean_bias_agg - mean_bias_day, 
                       breaks = c(-0.01, 0, 0.01, 0.02, 0.05),
                       labels = c("-0.01 - 0", "0 - 0.01", "0.01 - 0.02", 
                                  "> 0.02"))
  # Sharpness
  cat_diff_shar <- cut(mean_shar_agg - mean_shar_day, 
                       breaks = c(-0.05, -0.01, 0, 0.01, 1),
                       labels = c("-0.05 - -0.01", "-0.01 - 0", 
                                  "0 - 0.01", "> 0.01"))
  # RPS
  cat_diff_rps <- cut(mean_rps_agg - mean_rps_day, 
                      breaks = c(0, 0.005, 0.01, 0.02, 0.05),
                      labels = c("0 - 0.005", "0.005 - 0.01", 
                                 "0.01-0.02", "> 0.02"))
  ## Plot the three barplots representing the categories of difference
  par(mfrow = c(1,3), mar = c(5,5,4,1), las = 1, bty = "l", cex.axis = 1.1, 
      cex.main = 1.5, cex.lab = 2)
  barplot(table(cat_diff_bias), main = "", xlab = "Difference", 
          col = c(transp("grey", .6), transp("grey", .6),
                  transp("blue", .3), transp("blue", .6)),
          border = NA, ylab = "Count", ylim = c(0, 80))
  legend("top", fill = c(transp("purple", .9), transp("grey", .6),
                         transp("blue", .9)), border = NA,
         legend = c("Better in the aggregated model",
                    "Close in both models", 
                    "Better in the daily model"), bty = "n", cex = 1.5)
  mtext("A: Bias", side = 3, line = 1, adj = .1, cex = 2)
  barplot(table(cat_diff_shar), main = "", xlab = "Difference", 
          col = c(transp("purple", .3), transp("grey", .6),
                  transp("grey", .6), transp("blue", .3)),
          border = NA)
  mtext("B: Sharpness", side = 3, line = 1, adj = .1, cex = 2)
  barplot(table(cat_diff_rps), main = "", xlab = "Difference", 
          col = c(transp("grey", .3), transp("blue", .3),
                  transp("blue", .6), transp("blue", .9)),
          border = NA)
  mtext("C: RPS", side = 3, line = 1, adj = .1, cex = 2)
}
