figure_1 <- function(data, mean, sd, frac_no_miss = .5){
  ## Generate composite serial interval from mean, standard deviation and 
  ## proportion of connection without missing generation
  # Missing ancestor
  w0 <- dnorm(x = 1:50, mean = 0, sd = sqrt(2 * sd**2))
  # Direct connection
  w <- dnorm(x = 1:50, mean = mean, sd = sd)
  # Missing link
  w_dens2 <- conv(w, max_days = 50)
  # Composite SI
  w_composite <- w0 * (1-frac_no_miss)/2 + w*frac_no_miss + 
    w_dens2 * (1-frac_no_miss)/2
  
  # Intialise n_cases, the transmission potential for each day
  n_cases <- data * 0
  ## Compute n_cases throughout the time span
  t_min <- min(as.Date(rownames(data)))
  for (t in seq_len(nrow(data))[-1]){
    ## Convolve the number of cases in the last 50 days with the composite SI
    # Compute the dates 
    if(t<=50) t_i_min <- t_min else 
      t_i_min <- as.Date(rownames(data)[t - 50])
    date_t <- as.Date(rownames(data)[t])
    dates_t <- as.Date(t_i_min:(date_t-1),  origin = "1970-01-01")
    # Multiply the number of cases at each date with the reverse SI
    cases_tot <- data[as.character(dates_t),] * 
      matrix(w_dens[length(dates_t):1], nrow = length(dates_t), 
             ncol = ncol(data), byrow = F)
    # Compute n_cases at t
    n_cases[t-1,] <- colSums(cases_tot)
  }
  
  # Define the layout of the three plots (Two columns on the first row, 
  # one column on the second)
  layout_matrix <- matrix(1:4, ncol = 2, byrow = T)
  layout_matrix[2,2] <- 3
  layout(layout_matrix)
  par(mar = c(5, 5, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 2)
  # Plot the daily number of cases in the data
  plot(rowSums(data) ~ as.Date(rownames(data)), 
       type = "h", lend = "butt", ylab = "Number of cases", xlab = "Date")
  mtext("A", side = 3, line = -2, adj = 0.1, cex = 2)
  # Plot the composite serial interval
  plot(w_composite~seq_along(w_composite), type = "l", xlim = c(0, 30), 
       ylim = c(0, 0.15), xlab = "Days", ylab = "")
  title(ylab = "Probability", line = 3.5)
  mtext("B", side = 3, line = -2, adj = 0.1, cex = 2)
  # Add colours corresponding to the three "sections" of the composite SI
  lines(.25 * w0[w0>.025]~seq_along(w0[w0>.025]), type = "l", col = "orange",
        lwd = 2)
  lines(.5 * w[w>.025]~seq_along(w)[w>.025], type = "l", col = "red", lwd = 2)
  lines(.25 * w_dens2[w_dens2>.025]~seq_along(w_dens2)[w_dens2>.025], type = "l", 
        col = "darkred", lwd = 2)
  # Plot the transmission potential
  plot(rowSums(n_cases) ~ as.Date(rownames(data)), 
       type = "h", lend = "butt", ylab = "Transmission Potential", xlab = "Date")
  mtext("C", side = 3, line = -2, adj = 0.05, cex = 2)
  
}

figure_2 <- function(tot, nei, ar_labs = NULL, ne_labs = NULL, 
                     end_labs = NULL, other_labs = NULL){
  # Extract and sort the coefficients by components in the tot model
  ar_coef1 <- grepl("ar\\.", names(coef(tot))) %>% which
  ne_coef1 <- grepl("ne\\.", names(coef(tot))) %>% which
  end_coef1 <- grepl("end\\.", names(coef(tot))) %>% which
  other1 <- which(!is.element(seq_along(coef(tot)), 
                              c(ar_coef1, ne_coef1, end_coef1)))
  # Exclude the intercepts
  ar_coef1 <- ar_coef1[-1]
  ne_coef1 <- ne_coef1[-1]
  end_coef1 <- end_coef1[-1]
  
  # Extract and sort the coefficients by components in the neighbour model
  ar_coef2 <- grepl("ar\\.", names(coef(nei))) %>% which
  ne_coef2 <- grepl("ne\\.", names(coef(nei))) %>% which
  end_coef2 <- grepl("end\\.", names(coef(nei))) %>% which
  other2 <- which(!is.element(seq_along(coef(nei)), 
                              c(ar_coef2, ne_coef2, end_coef2)))
  
  # Exclude the intercepts
  ar_coef2 <- ar_coef2[-1]
  ne_coef2 <- ne_coef2[-1]
  end_coef2 <- end_coef2[-1]
  
  ## Plot the coefficients from the autoregressive component
  figs <- cbind(c(0, 0.65, 0, 0.5), c(0.5, 1, 0.5, 1),
                c(0.65, 1, 0, 0.5), c(0, 0.5, 0.5, 1))
  par(fig = figs[,4], mar = c(3, 5, 2, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 1.5)
  # from tot model
  y1 <- seq_along(ar_coef1) - .1
  plot(coef(tot)[ar_coef1] ~ y1, pch = 16, main = "", 
       ylab = "Effect", xaxt = "n", col = "blue", ylim = c(-.5, .5),
       xlim = c(y1[1] - .5, rev(y1)[1] + .7))
  # Add the 95% Confidence Intervals
  arrows(x0 = y1, y0 = confint(tot)[ar_coef1, 1],
         x1 = y1, y1 = confint(tot)[ar_coef1, 2], 
         angle = 90, code = 3, length = 0.05, col = "blue")
  mtext("A: Autoregressive component", side = 3, line = 0.5, adj = 0, cex = 1.5)
  # from nei model
  y2 <- seq_along(ar_coef2) + .1
  points(coef(nei)[ar_coef2] ~ y2, pch = 16, col = "purple")
  arrows(x0 = y2, y0 = confint(nei)[ar_coef2, 1],
         x1 = y2, y1 = confint(nei)[ar_coef2, 2], 
         angle = 90, code = 3, length = 0.05, col = "purple")
  # Add the legend
  legend("top", legend = c("Model 1: Spatial spread between all regions", 
                           "Model 2: Spatial spread between neighbours only"),
         lwd = 2, col = c("blue", "purple"), bty = "n", cex  = 1.1)
  # Add the axis labels
  if(!is.null(ar_labs)){
    axis(side = 1, at = seq_along(ar_coef1), labels = ar_labs)    
  } else
    axis(side = 1, at = seq_along(ar_coef1), 
         labels = gsub("ar.", "", names(coef(tot))[ar_coef1]))
  abline(h = 0, lty = 2)
  
  ## Plot the coefficients from the neighbourhood component
  # from tot model
  y1 <- seq_along(ne_coef1)
  if(!is.null(nei)) y1 <- y1 - .1
  par(fig = figs[,2], mar = c(3, 5, 2, 1), las = 1, bty ="l", new = T, cex.axis = 1.2,
      cex.main = 1.5, cex.lab = 1.5)
  plot(coef(tot)[ne_coef1] ~ y1, pch = 16, main = "", 
       ylab = "Effect", xaxt = "n", ylim = c(-.5, 1.2), col = "blue", 
       xlim = c(y1[1] - .5, rev(y1)[1] + .7))
  # Add the 95% Confidence Intervals
  arrows(x0 = y1, y0 = confint(tot)[ne_coef1, 1],
         x1 = y1, y1 = confint(tot)[ne_coef1, 2], 
         angle = 90, code = 3, length = 0.05, col = "blue")
  mtext("B: Neighbourhood component", side = 3, line = 0.5, adj = 0, cex = 1.5)
  
  # from nei model
  y2 <- seq_along(ne_coef2) + .1
  points(coef(nei)[ne_coef2] ~ y2, pch = 16, col = "purple")
  # Add the 95% Confidence Intervals
  arrows(x0 = y2, y0 = confint(nei)[ne_coef2, 1],
         x1 = y2, y1 = confint(nei)[ne_coef2, 2], 
         angle = 90, code = 3, length = 0.05, col = "purple")
  
  # Add the axis labels
  if(!is.null(ne_labs)){
    axis(side = 1, at = seq_along(ne_coef1), labels = ne_labs)
  } else
    axis(side = 1, at = seq_along(ne_coef1), 
         labels = gsub("ne.", "", names(coef(tot))[ne_coef1]))
  abline(h = 0, lty = 2)
  
  ## Plot the coefficients from the endemic component
  # from tot model
  y1 <- seq_along(end_coef1) - .1
  par(fig = figs[,1], mar = c(3, 5, 2, 1), las = 1, bty ="l", new = T, cex.axis = 1.2,
      cex.main = 1.5, cex.lab = 1.5)
  plot(coef(tot)[end_coef1] ~ y1, pch = 16, main = "", ylab = "Effect", 
       xaxt = "n", xlim = c(y1[1] - .5, rev(y1)[1] + .7), ylim = c(-1, 2), 
       col = "blue")
  mtext("C: Endemic component", side = 3, line = 0.5, adj = 0, cex = 1.5)
  # Add the 95% Confidence Intervals
  arrows(x0 = y1, y0 = confint(tot)[end_coef1, 1],
         x1 = y1, y1 = confint(tot)[end_coef1, 2], 
         angle = 90, code = 3, length = 0.05, col = "blue")
  
  # from nei model
  y2 <- seq_along(end_coef2) + .1
  points(coef(nei)[end_coef2] ~ y2, pch = 16, col = "purple")
  # Add the 95% Confidence Intervals
  arrows(x0 = y2, y0 = confint(nei)[end_coef2, 1],
         x1 = y2, y1 = confint(nei)[end_coef2, 2], 
         angle = 90, code = 3, length = 0.05, col = "purple")
  
  # Add the axis labels
  if(!is.null(end_labs)){
    axis(side = 1, at = seq_along(end_coef1), labels = end_labs)    
  } else
    axis(side = 1, at = seq_along(end_coef1), 
         labels = gsub("end.", "", names(coef(tot))[end_coef1]))
  abline(h = 0, lty = 2)
  
  ## Plot the other coefficients (gravity model and overdispersion)
  # from tot model
  y1 <- seq_along(other1) - .1
  par(fig = figs[,3], mar = c(3, 5, 2, 1), las = 1, bty ="l", new = T, cex.axis = 1.2,
      cex.main = 1.5, cex.lab = 1.5)
  plot(coef(tot)[other1] ~ y1, pch = 16, main = "", ylab = "Effect", 
       xaxt = "n", ylim = c(0, .8), col = "blue", 
       xlim = c(y1[1] - .5, rev(y1)[1] + .7))
  mtext("D: Others", side = 3, line = .5, adj = 0, cex = 1.5)
  # Add the 95% Confidence Intervals
  arrows(x0 = y1, y0 = confint(tot)[other1, 1],
         x1 = y1, y1 = confint(tot)[other1, 2], 
         angle = 90, code = 3, length = 0.05, col = "blue")
  # from nei model
  y2 <- seq_along(other2) + .1 + 1
  points(coef(nei)[other2] ~ y2, pch = 16, col = "purple")
  # Add the 95% Confidence Intervals
  arrows(x0 = y2, y0 = confint(nei)[other2, 1],
         x1 = y2, y1 = confint(nei)[other2, 2], 
         angle = 90, code = 3, length = 0.05, col = "purple")
  # Add the axis labels
  if(!is.null(other_labs)){
    axis(side = 1, at = seq_along(other1), labels = other_labs)    
  } else
    axis(side = 1, at = seq_along(other1), 
         labels = gsub("neweights.", "", names(coef(tot))[other1]))
  
    
}

figure_3 <- function(model, map){
  ### Compute the last value of the predictors
  ## Remove the seasonal effect (we are interested in the average value)
  # Compute the terms of the model
  terms_model <- terms(model)
  # Find out what columns of the term matrix generate the seasonal effect
  terms_season <- c(grep("cos", as.character(unlist(terms_model$terms[2,]))),
                    grep("sin", as.character(unlist(terms_model$terms[2,]))))
  # Set the effect of these columns to 0
  for (j in terms_season) 
    terms_model$terms[[1, j]] <- terms_model$terms[[1, j]] * 0
  # Compute the values of the predictors at each iteration
  mean_model <- meanHHH(coef(model), terms_model)
  ### Values of the endemic predictor
  map$type <- "Endemic"
  # Extract the value of the endemic predictor from mean_model
  map$Reff <- mean_model$end.exppred[nrow(mean_model$end.exppred),
                                     ][as.character(map$nuts3)]
  # Cut in categories
  map$R_cat <- cut(map$Reff, breaks = c(-1, 5e-4, 1e-3, 5e-3,
                                        max(map$Reff + .1, na.rm = T)), 
                   labels = c("<5e-4", "5-10e-4", "10-50e-4", ">5e-3"))
  map$lab <- ""
  # Generate the ggplot of the endemic predictor
  imp_end <- ggplot(map) +  geom_sf(aes(fill = R_cat)) + facet_grid(lab~type) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#feedde", "#fdae6b", "#e6550d", "#a63603"),
                      name = "Predictor") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), axis.ticks = element_blank(), 
          axis.line = element_blank(), legend.position = "bottom", 
          strip.text = element_text(size = 20))

  ## Values of the neighbourhood predictor
  map$type <- "Neighbourhood"
  # Extract the value of the neighbourhood predictor from mean_model
  map$Reff <- mean_model$ne.exppred[nrow(mean_model$ne.exppred),][
    as.character(map$nuts3)]
  # Normalise the values
  map$Reff <- map$Reff / 640
  # Cut in categories
  map$R_cat <- cut(map$Reff, breaks = c(-1, 0.25, 0.5, 1,
                                        max(map$Reff + .1, na.rm = T)), 
                   labels = c("0-0.25", "0.25-0.5", "0.5-1", ">1"))
  map$lab <- ""
  # Generate the ggplot of the neighbourhood predictor
  imp <- ggplot(map) +  geom_sf(aes(fill = R_cat)) + facet_grid(lab~type) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#feedde", "#fdae6b", "#e6550d", "#a63603"),
                      name = "Predictor") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), axis.ticks = element_blank(), 
          axis.line = element_blank(), legend.position = "bottom", 
          strip.text = element_text(size = 20))


  ## Values of the autoregressive predictor
  map$type <- "Autoregressive"
  # Extract the value of the autoregressive predictor from mean_model
  map$Reff <- mean_model$ar.exppred[nrow(mean_model$ar.exppred),][map$nuts3] + 
    mean_model$end.exppred[nrow(mean_model$ne.exppred),][map$nuts3]
  map$type <- "Autoregressive"
  break_vec <- c(-1, 0.65, 0.7, 0.75, max(map$Reff + .1, na.rm = T))
  break_lab <- c("0-0.65", "0.65-0.7", "0.7-0.75", ">0.75")
  
  # Cut in categories
  map$R_cat <- cut(map$Reff, breaks = break_vec, 
                   labels = break_lab)
  # Generate the ggplot of the autoregressive predictor
  in_reg_nb <- ggplot(map) +  geom_sf(aes(fill = R_cat)) + facet_grid(lab~type) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#feedde", "#fdae6b", "#e6550d", "#a63603"),
                      name = "Predictor") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), axis.ticks = element_blank(), 
          axis.line = element_blank(), legend.position = "bottom", 
          strip.text = element_text(size = 20))
  # Arrange all three plots together
  test_arrange <- arrangeGrob(imp_end, imp, in_reg_nb, ncol = 3)
  plot(test_arrange)
}

figure_4 <- function(){
  
}

figure_5_6 <- function(list_sim, model, map, which_dates, breaks,
                       labs, name_elements, thresh_cases, order_reg = NULL,
                       import_loc = NULL){
  ## For each figure, compute the proportion of simulation where the number
  ## of cases in each region was above the thresholds defined in thresh_cases
  list_maps <- lapply(list_sim, function(X){
    # Compute the number of cases per region in each simulation
    n_per_region_X <- matrix(lapply(X, function(Y)
      return(colSums(Y[which_dates,]))) %>% 
        unlist, nrow = length(X), byrow = T)
    colnames(n_per_region_X) <- colnames(model$stsObj@observed)
    # Initialise the three maps (per element in thresh_cases)
    map_X1 <- map_X2 <- map_X3 <- map
    ## For each region, compute the proportion of simulation where the number 
    ## of cases was above:
    # -thresh_cases[1]
    map_X1$sup_thresh <- 100 * 
      colSums(n_per_region_X >= thresh_cases[1])[as.character(map_X1$nuts3)]/
      nrow(n_per_region_X)
    # -thresh_cases[2]
    map_X2$sup_thresh <- 100 * 
      colSums(n_per_region_X >= thresh_cases[2])[as.character(map_X2$nuts3)]/
      nrow(n_per_region_X)
    # -thresh_cases[3]
    map_X3$sup_thresh <- 100 * 
      colSums(n_per_region_X >= thresh_cases[3])[as.character(map_X3$nuts3)]/
      nrow(n_per_region_X)
    # Turn the proportion of simulations into categories, defined in "breaks"
    map_X1$cat_thresh <- cut(map_X1$sup_thresh, breaks = breaks, labels = labs)
    map_X2$cat_thresh <- cut(map_X2$sup_thresh, breaks = breaks, labels = labs)
    map_X3$cat_thresh <- cut(map_X3$sup_thresh, breaks = breaks, labels = labs)
    # Create the label of each map
    if(thresh_cases[1] == 1){
      map_X1$cat <- paste0("At least ", thresh_cases[1], " case")
    } else map_X1$cat <- paste0("At least ", thresh_cases[1], " cases")
    map_X2$cat <- paste0("At least ", thresh_cases[2], " cases")
    map_X3$cat <- paste0("At least ", thresh_cases[3], " cases")
    return(rbind(map_X1, map_X2, map_X3))
  })
  # Name each element of list_maps with the names set in name_elements
  for (i in seq_len(length(list_maps))) list_maps[[i]]$type <- name_elements[i]
  
  ## Create the maps
  # Bind all the maps created in the same data frame
  map_tot <- do.call(rbind, list_maps)
  # If a particular order was given, format the column "type" as a factor
  if(!is.null(order_reg))
    map_tot$type <- factor(map_tot$type, levels = name_elements[order_reg])
  # Generate plot
  p_cases <- ggplot(map_tot) +  geom_sf(aes(fill = cat_thresh)) + 
    facet_grid(type~cat) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#2c7bb6", "#abd9e9", "#e0f3f8", 
                                 "#ffffbf", "#fdae61", "#d7191c"),
                      name = "Percentage of simulations") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + 
    theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), axis.ticks = element_blank(),  
          strip.text = element_text(size = 12),  axis.line = element_blank(), 
          legend.position = "bottom") + guides(fill = guide_legend(byrow = T))
  
  if(!is.null(import_loc)){
    map_tot$long_imp <- import_loc[as.character(map_tot$type),long]
    map_tot$lat_imp <- import_loc[as.character(map_tot$type),lat]
    p_cases <- p_cases + 
      geom_point(data = map_tot, aes(x = long_imp, y = lat_imp), size = 2) 
  }
  
  return(p_cases)
}

