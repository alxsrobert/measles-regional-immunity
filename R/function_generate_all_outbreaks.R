generate_all_sim <- function(n_sim, w_dens, cov, area, pop, list_params,
                             distance_matrix, start = as.Date("2009-01-29"),
                             rd_params = T){
  # Create a list where each element will be a matrix containing the number
  # of cases per region per day in a given simulation.
  list_sim <- vector(mode = "list", length = n_sim)
  # Compute the number of parameters
  n_params <- length(list_params$mean)
  # Update the names of the parameters
  names_params <- c("intercept_ar", "vacc_ar", "cat1_imm_ar", "cat2_imm_ar",
                    "pop_ar", "area_ar", "sin_ar", "cos_ar", 
                    "intercept_ne", "vacc_ne", "cat1_imm_ne", "cat2_imm_ne",
                    "pop_ne", "sin_ne", "cos_ne", 
                    "intercept_en", "vacc_en", "cat1_imm_en", "cat2_imm_en",
                    "pop_en", "sin_en", "cos_en", 
                    "dist", "pop", "overdisp")
  # Data frame of the parameters for each simulation
  params_sim <- as.data.frame(matrix(nrow = n_sim, ncol = n_params))
  # Initialise vector of parameters
  params_current <- numeric(n_params)
  colnames(params_sim) <- names(params_current) <- names_params
  # Sort the different data sources so that they have the same column order
  area <- area[, colnames(cov)]
  pop <- pop[, colnames(cov)]
  distance_matrix <- distance_matrix[colnames(cov), colnames(cov)]
  # 3 Data frames of the final value of the predictor for each region 
  r_fin_ar <- r_fin_ne <- r_fin_end <- 
    as.data.frame(matrix(nrow = n_sim, ncol = nrow(distance_matrix)))
  colnames(r_fin_ar) <- colnames(r_fin_ne) <- colnames(r_fin_end) <- 
    colnames(distance_matrix)
  # Matrix showing th proportion of cases coming from each component
  prop_comp <- matrix(nrow = n_sim, ncol = 3)
  
  for(i in 1:n_sim){
    # Draw the parameter sets using the covariance matrix
    if(rd_params == F) {
      dt_params <- matrix(list_params[[1]])
    } else {
      dt_params <- list_params[[1]] + (list_params[[2]] %>% chol %>% t)%*%
        matrix(rnorm(n = ncol(list_params[[2]])),
               nrow = ncol(list_params[[2]]), ncol = 1)
    }
    # Update params_current
    params_current[] <- dt_params[1:n_params,1]
    params_current["overdisp"] <- exp(params_current["overdisp"])
    # Run the simulation
    outbreak <- generate_sim(params = params_current, pop = pop, cov = cov, 
                             area = area, w_dens = w_dens, 
                             start = as.Date("2009-01-01"),  
                             distance_matrix = distance_matrix)
    # If more than 100,000 cases in ten years, re-simulate
    while(sum(outbreak$n_cases > 100000)){
      dt_params <- list_params[[1]] + (list_params[[2]] %>% chol %>% t)%*%
        matrix(rnorm(n = ncol(list_params[[2]])),
               nrow = ncol(list_params[[2]]), ncol = 1)
      params_current[] <- dt_params[1:n_params,1]
      
      params_current["overdisp"] <- exp(params_current["overdisp"])
      outbreak <- generate_sim(params = params_current, pop = pop, cov = cov, 
                               area = area, w_dens = w_dens, 
                               start = start,  
                               distance_matrix = distance_matrix)
    }
    # Extract the number of cases per day per region
    list_sim[[i]] <- outbreak$n_cases
    # Extract the final value of the autoregressive predictor per region
    r_fin_ar[i,] <- outbreak$r0_ar[nrow(outbreak$r0_ar),]
    # Extract the final value of the neighbourhood predictor per region
    r_fin_ne[i,] <- outbreak$r0_ne[nrow(outbreak$r0_ne),]
    # Extract the final value of the endemic predictor per region
    r_fin_end[i,] <- outbreak$r0_en[nrow(outbreak$r0_en),]
    # Compute the proportion of cases coming from each component
    all_cases <- c(AR = sum(outbreak$n_cases_ar, na.rm = T),
                   NE = sum(outbreak$n_cases_ne, na.rm = T),
                   EN = sum(outbreak$n_cases_en))
    prop_comp[i,] <- all_cases / sum(all_cases)
    # Extract the parameter set for this simulation
    params_sim[i,] <- params_current
  }
  
  return(list(sim = list_sim, r_ar = r_fin_ar, r_ne = r_fin_ne,
              r_end = r_fin_end, params_sim = params_sim,
              prop_comp = prop_comp))
}
