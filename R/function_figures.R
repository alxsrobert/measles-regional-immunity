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
                          exclude_overdisp = T){
  # If which_plot is not numeric, match the elements to the parameter names
  if(!is.numeric(which_plot)){
    which_plot <- which(is.element(colnames(list_all_sim$params_sim), which_plot))
  }
  # Generate CI matrices
  CI_day <- in_CI(list_all_sim = list_all_sim, hhh4_runs = hhh4_day, CI = CI)
  CI_agg <- in_CI(list_all_sim = list_all_sim, hhh4_runs = hhh4_agg, CI = CI)
  # Remove overdisp parameter
  if(exclude_overdisp){
    CI_day <- CI_day[-25,]
    CI_agg <- CI_agg[-25,]
  }
  # Generate difference between mean of the estimates and the parameter set
  dist_day <- distance_params(list_all_sim = list_all_sim, hhh4_runs = hhh4_day)
  dist_agg <- distance_params(list_all_sim = list_all_sim, hhh4_runs = hhh4_agg)
  # Generate pit values
  pit_i_day <- pit_param(list_all_sim = list_all_sim, hhh4_runs = hhh4_day)
  pit_i_agg <- pit_param(list_all_sim = list_all_sim, hhh4_runs = hhh4_agg)
  # Histograms of the proportion of parameters within the confidence interval
  par(mfrow = c(2,1), mar = c(3,4,2,0), bty = "l")
  hist(rowSums(CI_day)/ncol(CI_day), main = "", breaks = seq(0, 1, 0.025), 
       xlab = "")
  abline(v = .95, col = "red", lwd = 2, lty = 2)
  hist(rowSums(CI_agg)/ncol(CI_agg), main = "", breaks = seq(0, 1, 0.025), 
       xlab = "")
  abline(v = .95, col = "red", lwd = 2, lty = 2)
  
  ## For each parameter listed in which_plot, generate 3 figures:
  # 1- The density plot of the difference between estimates and parameter 
  # 2- Density plots of the mean of the parameter (with data)
  # 3- pit histogram of the figure
  par(mfrow = c(min(c(3, length(which_plot))),3), mar = c(5, 5, 1, 1), las = 1, 
      bty ="l", cex.axis = 1.1, cex.main = 1.5, cex.lab = 1.5)
  for(i in seq_along(which_plot)){
    param_i <- which_plot[i]
    # Generate density objects
    dens_dis_day <- density(dist_day[param_i,], na.rm = T)
    dens_dis_agg <- density(dist_agg[param_i,], na.rm = T)
    # Initialise empty plot
    plot(dens_dis_agg, col = "white", xlab = "", ylab = "", 
         main = colnames(list_all_sim$params_sim)[param_i], 
         ylim = c(0, max(c(dens_dis_agg$y, dens_dis_day$y))),
         xlim = c(max(-3, min(c(dens_dis_agg$x, dens_dis_day$x))),
                  max(c(dens_dis_agg$x, dens_dis_day$x))))
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
    plot(dens_agg, col = "white", xlab = "", ylab = "", 
         main = colnames(list_all_sim$params_sim)[param_i], 
         ylim = c(0, max(c(dens_day$y, dens_data$y, dens_agg$y))),
         xlim = c(max(c(min(c(dens_agg$x, dens_day$x, dens_data$x)), -3)), 
                  max(c(dens_agg$x, dens_day$x, dens_data$x))))
    # Add density plots
    polygon(dens_data, col=transp("black", .3), border=NA)
    polygon(dens_day, col=transp("#fc8d59", .3), border=NA)
    polygon(dens_agg, col=transp("#91bfdb", .3), border=NA)
    abline(v = 0, lty = 2, lwd = 2)
    # Generate histograms
    hist_day <- hist(pit_i_day[param_i,], breaks = seq(0,1, .1), plot = F)
    hist_agg <- hist(pit_i_agg[param_i,], breaks = seq(0,1, .1), plot = F)
    hist_day$counts <- hist_day$density
    hist_agg$counts <- hist_agg$density
    plot(hist_day, col = transp("#fc8d59", .3), ylab = "", border = NA,
         main = colnames(list_all_sim$params_sim)[param_i], 
         xlab = "", ylim = c(0, max(c(hist_day$counts, hist_agg$counts))))
    plot(hist_agg, col = transp("#91bfdb", .3), add = TRUE, border = NA)
    abline(h = 1, lwd = 2, lty = 2)
  }
}
