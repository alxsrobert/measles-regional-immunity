# Figure 1: Data description
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

# Figure 2: Parameter estimates
figure_2 <- function(tot, nei, ar_labs = NULL, ne_labs = NULL, 
                     end_labs = NULL, other_labs = NULL,
                     legend_text = c("Model 1: Spatial spread between all regions", 
                                     "Model 2: Spatial spread between neighbours only")){
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
  
  # Match the coefficients between the two models in each component
  pos_ar2 <- match(names(coef(nei))[ar_coef2], names(coef(tot))[ar_coef1])
  pos_ne2 <- match(names(coef(nei))[ne_coef2], names(coef(tot))[ne_coef1])
  pos_end2 <- match(names(coef(nei))[end_coef2], names(coef(tot))[end_coef1])
  pos_other2 <- match(names(coef(nei))[other2], names(coef(tot))[other1])
  
  if(any(is.na(pos_other2))) pos_other2[is.na(pos_other2)] <- 2
  
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
  y2 <- pos_ar2 + .1
  points(coef(nei)[ar_coef2] ~ y2, pch = 16, col = "purple")
  arrows(x0 = y2, y0 = confint(nei)[ar_coef2, 1],
         x1 = y2, y1 = confint(nei)[ar_coef2, 2], 
         angle = 90, code = 3, length = 0.05, col = "purple")
  # Add the legend
  legend("top", legend = legend_text, lwd = 2, col = c("blue", "purple"), 
         bty = "n", cex  = 1.1)
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
  y2 <- pos_ne2 + .1
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
  y2 <- pos_end2 + .1
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
  y2 <- pos_other2 + .1
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

# Figure 3: Geographical distribution of the predictors
figure_3 <- function(model1, model2, map){
  ### Compute the last value of the predictors
  ## Remove the seasonal effect (we are interested in the average value)
  # Compute the terms of each model
  terms_model1 <- terms(model1)
  terms_model2 <- terms(model2)
  
  # Find out what columns of the term matrix generate the seasonal effect
  terms_season1 <- c(grep("cos", as.character(unlist(terms_model1$terms[2,]))),
                     grep("sin", as.character(unlist(terms_model1$terms[2,]))))
  terms_season2 <- c(grep("cos", as.character(unlist(terms_model2$terms[2,]))),
                     grep("sin", as.character(unlist(terms_model2$terms[2,]))))
  # Set the effect of these columns to 0
  for (j in terms_season1) 
    terms_model1$terms[[1, j]] <- terms_model1$terms[[1, j]] * 0
  for (j in terms_season2) 
    terms_model2$terms[[1, j]] <- terms_model2$terms[[1, j]] * 0
  # Compute the values of the predictors at each iteration
  mean_model1 <- meanHHH(coef(model1), terms_model1)
  mean_model2 <- meanHHH(coef(model2), terms_model2)
  ### Values of the endemic predictor
  map$type <- "Endemic"
  map2 <- map
  map$model <- "Model 1"
  map2$model <- "Model 2"
  # Extract the value of the endemic predictor from mean_model
  map$Reff <- mean_model1$end.exppred[nrow(mean_model1$end.exppred),
  ][as.character(map$nuts3)]
  map2$Reff <- mean_model2$end.exppred[nrow(mean_model2$end.exppred),
  ][as.character(map2$nuts3)]
  # Cut in categories
  map$R_cat <- cut(map$Reff, breaks = c(-1, 5e-4, 1e-3, 5e-3,
                                        max(map$Reff + .1, na.rm = T)), 
                   labels = c("<5e-4", "5-10e-4", "10-50e-4", ">5e-3"))
  map2$R_cat <- cut(map2$Reff, breaks = c(-1, 5e-4, 1e-3, 5e-3,
                                          max(map2$Reff + .1, na.rm = T)), 
                    labels = c("<5e-4", "5-10e-4", "10-50e-4", ">5e-3"))
  maptot <- rbind(map, map2)
  # Generate the ggplot of the endemic predictor
  imp_end <- ggplot(maptot) +  geom_sf(aes(fill = R_cat)) + facet_grid(model~type) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#feedde", "#fdae6b", "#e6550d", "#a63603"),
                      name = "Predictor") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), axis.ticks = element_blank(), 
          axis.line = element_blank(), legend.position = "bottom", 
          strip.text = element_text(size = 20),
          strip.text.y = element_blank())
  
  ## Values of the neighbourhood predictor
  map$type <- "Neighbourhood"
  map2$type <- "Neighbourhood"
  # Extract the value of the neighbourhood predictor from mean_model
  map$Reff <- mean_model1$ne.exppred[nrow(mean_model1$ne.exppred),][
    as.character(map$nuts3)]
  map2$Reff <- mean_model2$ne.exppred[nrow(mean_model2$ne.exppred),][
    as.character(map2$nuts3)]
  # Normalise the values
  map$Reff <- map$Reff / 300
  map2$Reff <- map2$Reff / 300
  # Cut in categories
  map$R_cat <- cut(map$Reff, breaks = c(-1, 0.25, 0.5, 1,
                                        max(map$Reff + .1, na.rm = T)), 
                   labels = c("0-0.25", "0.25-0.5", "0.5-1", ">1"))
  map2$R_cat <- cut(map2$Reff, breaks = c(-1, 0.25, 0.5, 1,
                                          max(map2$Reff + .1, na.rm = T)), 
                    labels = c("0-0.25", "0.25-0.5", "0.5-1", ">1"))
  maptot <- rbind(map, map2)
  # Generate the ggplot of the neighbourhood predictor
  imp <- ggplot(maptot) +  geom_sf(aes(fill = R_cat)) + facet_grid(model~type) +
    scale_fill_manual(na.translate = F, guide = guide_legend(),
                      values = c("#feedde", "#fdae6b", "#e6550d", "#a63603"),
                      name = "Predictor") + 
    coord_sf(ylim = c(42, 51.8), xlim = c(-4.5, 7.7)) + theme_classic(base_size = 12) + 
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          strip.background = element_blank(), axis.ticks = element_blank(), 
          axis.line = element_blank(), legend.position = "bottom", 
          strip.text = element_text(size = 20),
          strip.text.y = element_blank())
  
  
  ## Values of the autoregressive predictor
  map$type <- "Autoregressive"
  map2$type <- "Autoregressive"
  # Extract the value of the autoregressive predictor from mean_model
  map$Reff <- mean_model1$ar.exppred[nrow(mean_model1$ar.exppred),][map$nuts3] + 
    mean_model1$end.exppred[nrow(mean_model1$ne.exppred),][map$nuts3]
  map2$Reff <- mean_model2$ar.exppred[nrow(mean_model2$ar.exppred),][map2$nuts3] + 
    mean_model2$end.exppred[nrow(mean_model2$ne.exppred),][map2$nuts3]
  map$type <- "Autoregressive"
  map2$type <- "Autoregressive"
  break_vec <- c(-1, 0.65, 0.7, 0.75, max(map$Reff + .1, na.rm = T))
  break_lab <- c("0-0.65", "0.65-0.7", "0.7-0.75", ">0.75")
  
  # Cut in categories
  map$R_cat <- cut(map$Reff, breaks = break_vec, 
                   labels = break_lab)
  map2$R_cat <- cut(map2$Reff, breaks = break_vec, 
                    labels = break_lab)
  maptot <- rbind(map, map2)
  # Generate the ggplot of the autoregressive predictor
  in_reg_nb <- ggplot(maptot) +  geom_sf(aes(fill = R_cat)) + facet_grid(model~type) +
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
  return(test_arrange)
}

# Figure 4: Fit and calibration
figure_4 <- function(model, list_calib){
  # Generate the number of cases stemming of each component at each date
  hhh_mod <- meanHHH(coef(model), terms(model))
  
  # Organise the plotting area
  layout_matrix <- matrix(rep(seq_len(4), each = 2), ncol = 4, byrow = T)
  layout_matrix[2,] <- 3:6
  layout(layout_matrix)
  par(mar = c(5,4,1,0), las = 1, bty = "l", cex.axis = 1.1, cex.main = 1.5, 
      cex.lab = 2)
  ## First plot: Daily fit
  # y = The observed data
  y <- model$control$data$response[-1,] %>% rowSums()
  dates <- (as.Date(names(y)))
  # Initialise an empty plot
  plot(y~dates, type = "n", xlab = "", ylab = "")
  # Plot the overall number of cases at each date
  lines((hhh_mod$epidemic + hhh_mod$endemic) %>% rowSums()~dates,
        col = "orange", type = "h")
  # Plot the autoregressive + endemic
  lines((hhh_mod$epi.own + hhh_mod$endemic) %>% rowSums()~dates,
        col = "blue", type = "h")
  # Plot the endemic
  lines(hhh_mod$endemic %>% rowSums()~dates, col = "grey", type = "h")
  # Add the data points
  points(y[y>0]~dates[y>0], type = "p", pch = 19, cex = .3, col = transp("black", .3))
  # Add the title
  title(xlab = "Days", ylab = "Number of cases", line = 2.5)
  legend("top", lwd = c(2,2,2,NA), bty = "n", cex = 1.5, 
         pch = c(NA, NA, NA, 19), col = c("orange", "blue", "grey", transp("black", .5)),
         legend = c("Neighbourhood component", "Autoregressive component",
                    "Endemic component", "Data"))
  mtext("A", side = 3, line = -2, adj = .1, cex = 2)
  
  ## Weekly fits
  par(mar = c(5, 2, 1, 2), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 2)
  # Generate the dates of weekly aggregation
  group_r <- seq(1, ceiling(nrow(hhh_mod$epi.own)/7))
  group_r <- rep(group_r, each = 7)[1:nrow(hhh_mod$epi.own)]
  
  
  # Aggregate the autoregressive, neighbourhood, and endemic components
  ar_agg <- aggregate.data.frame(x = hhh_mod$epi.own, by = list(group_r),
                                 FUN = sum)[,-1] %>% rowSums()
  ne_agg <- aggregate.data.frame(x = hhh_mod$epi.neighbours, by = list(group_r),
                                 FUN = sum)[, -1] %>% rowSums()
  end_agg <- aggregate.data.frame(x = hhh_mod$endemic, by = list(group_r),
                                  FUN = sum)[, -1] %>% rowSums()
  # Aggregate the reported data
  y <- aggregate.data.frame(x = model$control$data$response[-1,],
                            by = list(group_r), FUN = sum)[,-1] %>% rowSums()
  # Total aggregated number of cases
  tot_agg <- ar_agg + ne_agg + end_agg
  # Aggregated autoregressive + endemic
  ar_end_agg <- ar_agg + end_agg
  # Compute the unique dates of aggregation
  rownames(hhh_mod$epi.own) <- rownames(model$control$data$response[-1,])
  dates <- (as.Date(rownames(hhh_mod$epi.own[which(!duplicated(group_r)),])))
  ## Plot the aggregated values
  plot(y~dates, type = "n", xlab = "", ylab = "")
  lines((tot_agg)~dates, col = "orange", lwd = 3, type = "h")
  lines((ar_end_agg)~dates, col = "blue", lwd = 3, type = "h")
  lines(end_agg~dates, col = "grey", lwd = 3, type = "h")
  points(y[y>0]~dates[y>0], type = "p", pch = 19, cex = .3, col = transp("black", .5))
  title(xlab = "Weeks", line = 2.5)
  mtext("B", side = 3, line = -2, adj = 0.1, cex = 2)
  
  ### Generate the PIT histograms
  ## Define the function to compute the PIT values
  par(mar = c(5, 4, 1, 1), las = 1, bty ="l", cex.axis = 1.1,
      cex.main = 1.5, cex.lab = 2)
  Fbar1 <- function (u, Px, Pxm1){
    F_u <- punif(u, Pxm1, Px)  # also works for Pxm1 == Px => F_u = u >= Pxm1
    mean(F_u)
  }
  
  # For each element of list_calib, compute and plot the PIT histograms
  for(i in seq_along(list_calib)){
    breaks <- (0:10)/10
    calib <- list_calib[[i]]
    # Compute proportion of entry in each category of the PIT histogram
    Fbar_seq_prop <- vapply(X = breaks, FUN = Fbar1, FUN.VALUE = 0,
                            Px = c(calib$px),Pxm1 = c(calib$pxm1),
                            USE.NAMES = FALSE)
    # Plot the PIT histograms
    plot(10 * diff(Fbar_seq_prop)~breaks[-length(breaks)],
         type = "h", ylim = c(0.5, 1.5), lend = "butt", lwd = 20, col = "grey",
         ylab = "", xlab = "")
    abline(col = "red", lty = 2, h = 1)
    abline(col = "red", lty = 3, h = c(0.9, 1.1))
    title(ylab = "Relative Frequency", line = 2.5)
    # Add legends, titles, and labels
    if(i == 1){
      mtext("C: 3 days", side = 3, line = -2, adj = 0.1, cex = 2)
      par(mar = c(5, 2, 1, 2), las = 1, bty ="l", cex.axis = 1.1,
          cex.main = 1.5, cex.lab = 2)
    }
    if(i == 2) mtext("D: 7 days", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i == 3) mtext("E: 10 days", side = 3, line = -2, adj = 0.1, cex = 2)
    if(i == 4) mtext("F: 14 days", side = 3, line = -2, adj = 0.1, cex = 2)
  }
  title(xlab = "Probability Integral Transform", outer = T, line = - 1.5)
}

# Figure 5: 1-year-ahead simulations
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

