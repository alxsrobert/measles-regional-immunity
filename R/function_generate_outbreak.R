generate_sim <- function(w_dens, cov, area, pop, distance_matrix, 
                         params, start, thresh_incid = c(10, 45)){
  # Set up intercept of each component
  R0_ar <- rep(params["intercept_ar"], ncol(cov))
  R0_en <- rep(params["intercept_en"], ncol(cov))
  R0_ne <- rep(params["intercept_ne"], ncol(cov))
  names(R0_ar) <- names(R0_en) <- names(R0_ne) <- colnames(cov)
  
  ## Compute distance matrix
  if(!is.matrix(pop)){
    # If the number of inhabitants is constant throughout the study period
    # Exponential gravity model
    wji <- exp(-distance_matrix*(params["dist"])) * (pop^params["pop"]) / pop
    # Distance between a region and itself is 0
    wji[distance_matrix == 0] <- 0
  } else{
    # If the number of inhabitants varies throughout the study period
    # Exponential gravity model applied for each value of pop
    wji <- apply(pop, 1, function(X){
      if(any(names(params) == "dist")){
        wji_X <- exp(-distance_matrix*(params["dist"])) * (X^params["pop"]) / X 
      } else {
        wji_X <- (X^params["pop"]) / X
      }
      # Distance between a region and itself is 0
      wji_X[distance_matrix == 0] <- 0
      return(list(wji_X))})
    # REturn the first value to avoid having a list within a list
    wji <- lapply(wji, function(X) return(X[[1]]))
  }
  # Extract the first value for each covariate
  cov_0 <- cov[as.character(start), ]
  pop_0 <- pop[as.character(start),]
  area_0 <- area[as.character(start),]
  
  # At the starting point, the level of recent incidence is 0 in each region 
  incid_0 <- rep(1, length(cov_0))
  incid_1 <- rep(0, length(cov_0))
  incid_2 <- rep(0, length(cov_0))
  names(incid_0) <- names(incid_1) <- names(incid_2) <- names(cov_0)

  # Compute the number of iterations
  weeks_sim <- rownames(cov_ts)[rownames(cov_ts) >= start]
  # Convert to numeric
  num_dates <- as.numeric(as.Date(weeks_sim))
  
  # Initialise the matrix containing the number of cases per day per region
  n_cases <- matrix(nrow= length(weeks_sim), ncol = ncol(cov))
  rownames(n_cases) <- weeks_sim
  colnames(n_cases) <- colnames(cov)
  
  # Initialise the matrices with the values of the predictor per day per region
  all_r0_ar <- matrix(nrow = length(weeks_sim), ncol = ncol(cov))
  all_r0_ne <- matrix(nrow = length(weeks_sim), ncol = ncol(cov))
  all_r0_en <- matrix(nrow = length(weeks_sim), ncol = ncol(cov))
  
  # Compute the first value of each predictor
  t <- 0
  R0_ar_reg <- exp(R0_ar + params["vacc_ar"] * log(1 - cov_0/100) + 
                     params["cat1_imm_ar"] * incid_1 + 
                     params["cat2_imm_ar"] * incid_2 + 
                     params["pop_ar"] * log(pop_0 / 1000000) + 
                     params["area_ar"] * log(area_0) + 
                     params["sin_ar"] * sin(t * 2*pi/365)+
                     params["cos_ar"] * cos(t * 2*pi/365))
  R0_ne_reg <- exp(R0_ne + params["vacc_ne"] * log(1 - cov_0/100) + 
                     params["cat1_imm_ne"] * incid_1 + 
                     params["cat2_imm_ne"] * incid_2 + 
                     params["pop_ne"] * log(pop_0 / 1000000) + 
                     params["sin_ne"] * sin(t * 2*pi/365)+
                     params["cos_ne"] * cos(t * 2*pi/365))
  R0_en_reg <- exp(R0_en + params["vacc_en"] * log(1 - cov_0/100) + 
                     params["cat1_imm_en"] * incid_1 + 
                     params["cat2_imm_en"] * incid_2 + 
                     params["pop_en"] * log(pop_0 / 1000000) + 
                     params["sin_en"] * sin(t * 2*pi/365)+
                     params["cos_en"] * cos(t * 2*pi/365))
  # Add the initial values of the predictors to the predictor matrices
  all_r0_ar[1,] <- R0_ar_reg
  all_r0_ne[1,] <- R0_ne_reg
  all_r0_en[1,] <- R0_en_reg
  
  # Number of cases at t=0 is 0
  n_cases_ar <- n_cases_en <- n_cases_ne <- incid <- n_cases
  # Draw the number of cases at t = 1
  n_cases[t + 1,] <- sapply(R0_en_reg, function(X)
    rnbinom(1, mu = X, size = params["overdisp"]))
  # All cases at t=1 come from the endemic component
  n_cases_en[t + 1,] <- R0_en_reg
  
  for(i in seq_along(weeks_sim)[-1]){
    # Increase t by 1
    t <- t + 1
    time_i <- rownames(n_cases)[i]
    # Update the value of the covariates
    cov_i <- cov[as.character(time_i), ]
    area_i <- area[as.character(time_i),]
    pop_i <- pop[as.character(time_i),]

    ## Compute the level of recent transmission:
    if(as.Date(time_i) - as.Date(weeks_sim[1]) > 30){
      # Compute the numeric value of the current date since numeric comparison 
      # is faster
      num_date <- as.numeric(as.Date(time_i))
      ## For each region, sum cases infected 1 month to 1 year before time_i
      # generate dates
      weeks_incid <- 
        rownames(n_cases)[num_dates < num_date - 30 &
                            num_dates > num_date - 365 * 3]
      # Select rows
      incidence <- n_cases[weeks_incid,]
      # Sum rows
      if(length(weeks_incid) > 1) incidence <- colSums(n_cases[weeks_incid,]) 
      
      # Compute the number of cases per million
      incidence <- incidence / pop_i * 1000000
      incid_0 <- incid_1 <- incid_2 <- incidence * 0
      # Region with a number of cases per million below thresh_incid[1] belong to cat0
      incid_0[incidence < thresh_incid[1]] <- 1
      # Region with a number of cases per million thresh_incid[1]-thresh_incid[2] belong to cat1
      incid_1[incidence > thresh_incid[1] & incidence <= thresh_incid[2]] <- 1
      # Region with a number of cases per million above thresh[2] belong to cat2
      incid_2[incidence > thresh_incid[2]] <- 1
      if(any(incid_0 + incid_1 + incid_2 !=1)) stop("Problem incidence")
    }
    # Compute the value of the predictors at time_i
    R0_ar_reg <- exp(R0_ar +
                       params["vacc_ar"] * log(1 - cov_i/100) + 
                       params["cat1_imm_ar"] * incid_1 + 
                       params["cat2_imm_ar"] * incid_2 + 
                       params["pop_ar"] * log(pop_i / 1000000) + 
                       params["area_ar"] * log(area_i) + 
                       params["sin_ar"] * sin(t * 2*pi/365)+
                       params["cos_ar"] * cos(t * 2*pi/365))
    R0_ne_reg <- exp(R0_ne +
                       params["vacc_ne"] * log(1 - cov_i/100) + 
                       params["cat1_imm_ne"] * incid_1 + 
                       params["cat2_imm_ne"] * incid_2 + 
                       params["pop_ne"] * log(pop_i / 1000000) + 
                       params["sin_ne"] * sin(t * 2*pi/365)+
                       params["cos_ne"] * cos(t * 2*pi/365))
    R0_en_reg <- exp(R0_en +
                       params["vacc_en"] * log(1 - cov_i/100) + 
                       params["cat1_imm_en"] * incid_1 + 
                       params["cat2_imm_en"] * incid_2 + 
                       params["pop_en"] * log(pop_i/ 1000000) + 
                       params["sin_en"] * sin(t * 2*pi/365)+
                       params["cos_en"] * cos(t * 2*pi/365))
    # Update the matrices
    all_r0_ar[i,] <- R0_ar_reg
    all_r0_ne[i,] <- R0_ne_reg
    all_r0_en[i,] <- R0_en_reg
    ## Compute the infection potential using the serial interval and the number
    ## of cases in the last 50 days
    diff_time <- (seq_len(min(50, i - 1))) + max(0, i - 1 -50)
    cases_tot <- (n_cases[diff_time,] * 
                    matrix(w_dens[length(diff_time):1], 
                           nrow = length(diff_time), ncol = ncol(n_cases), 
                           byrow = F)) %>% colSums()
    names(cases_tot) <- colnames(n_cases)
    # Compute the average number of new cases per component
    new_ar_av <- R0_ar_reg * cases_tot
    new_en_av <- R0_en_reg
    if(!is.list(wji)) new_ne_av <- R0_ne_reg * ((cases_tot * wji) %>% colSums()) else
      new_ne_av <- R0_ne_reg * ((cases_tot * wji[[i-1]]) %>% colSums())
    # Draw the number of cases at time_i per region
    n_cases[i,] <- sapply(new_en_av + new_ar_av + new_ne_av, function(X)
      rnbinom(1, mu = X, size = params["overdisp"]))
    # Update the matrices of the average number of cases per component
    n_cases_ar[i,] <- new_ar_av
    n_cases_en[i,] <- new_en_av
    n_cases_ne[i,] <- new_ne_av
    # Update the matrix of the level of coverage per day per region
    incid[i,] <- incid_1 + incid_2 * 2
  }
  
  return(list(n_cases = n_cases, incid = incid,
              r0_ar = all_r0_ar, r0_ne = all_r0_ne, r0_en = all_r0_en, 
              n_cases_ar = n_cases_ar, n_cases_en = n_cases_en, 
              n_cases_ne = n_cases_ne))
}
