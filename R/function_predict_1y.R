## Simulate lists of outbreaks
simulation_main_figures <- function(model, n_param, n_sim_param, days, 
                                    data, pop_mat, w_dens, reg_import = NULL, 
                                    dates_import = NULL, thresh = NA){
  # Total number of simulations is equal to the number of parameter sets * 
  # the number of simulation per parameter set
  n_sim <- n_param * n_sim_param 
  # Extract the covariance matrix
  cov_mat <- model$cov
  
  # Compute and randomly draw the parameter sets
  dt_params <- model$coefficients + (cov_mat %>% chol %>% t)%*%
    matrix(rnorm(n = ncol(cov_mat)*n_sim),
           nrow = ncol(cov_mat), ncol = n_sim)
  # Initialise the list of simulations
  list_sim1 <- list_sim2 <- list_sim3 <- list_sim4 <- as.list(rep(NA, n_sim))
  
  ## Generate the different conditions of simulation:
  # 1: Increase the vaccination coverage in each region by 3%
  unvax_better <- log((exp(model$control$data$unvax)[nrow(model$control$data$unvax),] -
                         .03)) %>% as.matrix %>% t
  unvax_better[is.na(unvax_better) | unvax_better < log(.01)] <- log(.01)
  unvax_better <- matrix(unvax_better, nrow = length(days), ncol = ncol(unvax_better),
                         byrow = T)
  
  # 2: Decrease the vaccination coverage in each region by 3%
  unvax_worse <- log((exp(model$control$data$unvax)[nrow(model$control$data$unvax),] +
                        .03)) %>% as.matrix %>% t
  unvax_worse <- matrix(unvax_worse, nrow = length(days), ncol = ncol(unvax_worse),
                        byrow = T)
  # 3: Set the level of recent incidence to minimum in each region 
  elim <- (exp(model$control$data$cat1)[nrow(model$control$data$unvax),] * 0) %>% 
    as.matrix %>% t
  
  # Use the data up to the first day of prediction
  observed_sim <- data[rownames(data) < min(days),]
  
  for(j in seq_len(n_param)){
    # Initialise model_iter, and set the coefficients to the i-th parameter set
    model_iter <- model
    model_iter$coefficients <- dt_params[,j]
    # Generate one year of future outbreaks 
    set.seed(1)
    sim <- generate_1y_outbreak(model = model_iter, thresh = thresh, 
                                pop_mat = pop_mat, time = days, w_dens = w_dens,
                                n_sim = n_sim_param, bool_incidence = T,
                                data = observed_sim, reg_import = reg_import,
                                dates_import = dates_import)
    # Add these simulations to list_sim1
    list_sim1[seq_len(n_sim_param) + n_sim_param * (j-1)] <- sim
    
    # Generate one year of future outbreaks with increase of vaccine coverage
    set.seed(1)
    sim <- generate_1y_outbreak(model = model_iter, thresh = thresh,
                                pop_mat = pop_mat, time = days, w_dens = w_dens, 
                                data = observed_sim, n_sim = n_sim_param, 
                                bool_incidence = T, reg_import = reg_import,
                                conditions = list(unvax = unvax_better),
                                dates_import = dates_import)
    # Add these simulations to list_sim2
    list_sim2[seq_len(n_sim_param) + n_sim_param * (j-1)] <- sim
    
    # Generate one year of future outbreaks with decrease of vaccine coverage
    set.seed(1)
    sim <- generate_1y_outbreak(model = model_iter, thresh = thresh,
                                pop_mat = pop_mat, time = days, w_dens = w_dens, 
                                data = observed_sim, reg_import = reg_import, 
                                n_sim = n_sim_param, bool_incidence = T,
                                conditions = list(unvax = unvax_worse),
                                dates_import = dates_import)
    # Add these simulations to list_sim3
    list_sim3[seq_len(n_sim_param) + n_sim_param * (j-1)] <- sim
    
    # Generate one year of future outbreaks if no transmission in the past 3y
    set.seed(1)
    sim <- generate_1y_outbreak(model = model_iter, thresh = thresh,
                                pop_mat = pop_mat, time = days, w_dens = w_dens, 
                                data = observed_sim, reg_import = reg_import, 
                                n_sim = n_sim_param, bool_incidence = F,
                                dates_import = dates_import)
    # Add these simulations to list_sim4
    list_sim4[seq_len(n_sim_param) + n_sim_param * (j-1)] <- sim
    
  }
  list_sim <- list(list_sim1, list_sim2, list_sim3, list_sim4)
  return(list_sim)
}

## Simulate 1year of outbreak
generate_1y_outbreak <- function(model, pop_mat, time, w_dens, n_sim, 
                                 bool_incidence, data, conditions = list(),
                                 reg_import = NULL, dates_import = NULL, 
                                 thresh = NA){
  # Set the dates of the observations
  date_obs <- as.Date(rownames(data))
  # Set the dates of prediction
  date_ncas <- time
  
  ## Compute threshold for incidence category
  if(is.na(thresh)){
    
    data <- model$control$data$response
    all_dates <- as.Date(rownames(data))
    ## Compute category of incidence
    cases_3years <- t(sapply(data %>% rownames %>% as.Date, function(X){
      # Compute the number of cases reported 1 month to  3 years before X
      if(X < all_dates %>% min + 30){
        incidence <- data[1:2,]*0
      } else if(X < all_dates %>% min + 365 *3){
        incidence <- data[all_dates < X - 30,]
      } else{
        incidence <- data[all_dates < X - 30 &
                            all_dates > X - 365 * 3,]
      }
      if(!is.matrix(incidence)) res <- incidence else res <- colSums(incidence)
      return(res)
    }))
    rownames(cases_3years) <- data %>% rownames
    # Compute the number of cases per million
    incidence <- cases_3years / pop_mat[rownames(cases_3years), ] * 1000000
    
    thresh <- quantile(incidence, 2/3)
  }
  
  # Initialise the number of daily cases per region
  n_cases <- matrix(0, nrow = length(time), ncol = model$nUnit * n_sim)
  rownames(n_cases) <- date_ncas %>% as.character
  colnames(n_cases) <- rep(colnames(data), n_sim)
  # Format the composite serial interval as a matrix (for each region)
  mat_w <- matrix(rev(w_dens),nrow = length(w_dens), 
                  ncol = ncol(data), byrow = F)
  ## Compute the level of recent incidence
  number_recent <- colSums(data[date_obs > (time[1] - 365*3) &
                                  date_obs < (time[1] - 30),])
  # If scenario without recent incidence, set number_recent to 0
  if(bool_incidence == F) number_recent <- number_recent * 0
  ## Compute the category of recent incidence
  incidence_0 <- incidence_1 <- incidence_2 <- number_recent * 0
  incidence <- number_recent / pop_mat[nrow(pop_mat),] * 1000000
  incidence_0[incidence <= 10] <- 1
  incidence_1[incidence > 10 & incidence <= thresh] <- 1
  incidence_2[incidence > thresh] <- 1
  # Compute the connectivity matrix
  wji <- getNEweights(model)
  wji <- wji[,, dim(wji)[3]]
  rownames(wji) <- colnames(wji) <- colnames(data)
  # extract condition of simulation if any (eg changes in the vaccine coverage)
  conditions_0 <- lapply(conditions, function(X) return(X[1,]))
  ### Compute the values of predictors for the days of simulation using the 
  ### latest values of population, vaccine coverage, and surface 
  ### (excluding recent incidence).
  # Last date of observation
  max_t <- max(date_obs)
  # Row number of the last date of observation
  min_pred <- model$control$data$t[which(rownames(model$control$data$unvax) ==
                                           max_t)]
  # Dates of prediction
  dates_pred <- seq_along(time) + min_pred
  # Set the dates of prediction
  model$control$data$t <- dates_pred
  
  # Set all the covariates in data to their last value, with the number of rows
  # corresponding to the number of dates predicted
  model$control$data <- lapply(model$control$data, function(X){
    if(is.matrix(X)){
      X <- matrix(X[nrow(X),], nrow = length(date_ncas),
                  ncol = ncol(X), byrow = T)
      rownames(X) <- as.character(time)
    }
    return(X)
  })
  # Set the level of recent incidence to the minimum level 
  # (it will be included) in the calculation of the predictors later on, at 
  # this stage we only include covariates that do not depend on the nb of cases
  model$control$data$cat1 <- model$control$data$cat1 * 0
  model$control$data$cat2 <- model$control$data$cat2 * 0
  # Need to have the same number of rows in the elements of stsObj as in 
  # the elements of data for terms() to work. But it has no impact on the 
  # values of the predictors.
  model$stsObj <- model$stsObj[seq_along(dates_pred),]
  model$control$subset <- 2:(nrow(model$fitted.values)+1)
  
  # Compute the function terms with the new conditions
  terms_pred <- term_predict(model = model, conditions = conditions_0, 
                             t = dates_pred)
  
  # Compute the last value of transmission potential
  FI_t <- colSums(data[(nrow(data) - length(w_dens)+1) : nrow(data),] * mat_w)
  
  # Compute the value of the predictors (excluding the impact of incidence)
  mean_t <- meanHHH(model$coefficients, terms_pred)
  # Compute the values of the predictors adding in the impact of incidence
  r0_ar1 <- mean_t$ar.exppred[1,] * exp(incidence_1 * coef(model)["ar.cat1"]) * 
    exp(incidence_2 * coef(model)["ar.cat2"])
  r0_ne1 <- mean_t$ne.exppred[1,] * exp(incidence_1 * coef(model)["ne.cat1"]) * 
    exp(incidence_2 * coef(model)["ne.cat2"])
  r0_en1 <- mean_t$end.exppred[1,] * exp(incidence_1 * coef(model)["end.cat1"]) * 
    exp(incidence_2 * coef(model)["end.cat2"])
  # Compute the average number of new cases per compartment
  new_ar_av <- r0_ar1 * FI_t 
  new_en_av <- r0_en1 
  new_ne_av <- r0_ne1 * ((FI_t * wji) %>% colSums()) 
  if(!is.null(reg_import)){
    reg_import_nb <- which(is.element(colnames(n_cases), reg_import))
    new_ar_av <- 0
    new_ne_av <- 0
    new_en_av <- 0
    n_cases[1,][reg_import_nb] <- sum(dates_import == time[1])
  } else
    # Generate the number of new cases per region
    n_cases[1,] <- rnbinom(n_sim * length(new_ar_av), mu = new_en_av,
                           size = 1/(coef(model)["overdisp"]))
  ### Same thing for each time step
  for (i in seq_along(time)[-1]){
    # Compute the number of cases per region in the last 3 years (from data)
    incidence_i <- colSums(data[date_obs > (time[i] - 365*3) &
                                  date_obs < (time[i] - 30),])
    # If no recent incidence, set incidence_i to 0
    if(bool_incidence == F) incidence_i <- incidence_i * 0
    # Compute the number of new cases per region in the simulated time steps
    new_incidence <- n_cases[date_ncas > (time[i] - 365*3) &
                               date_ncas < (time[i] - 30),]
    # Add up incidence_i and new_incidence
    if(!is.matrix(new_incidence) && length(new_incidence) > 0){ 
      incidence_i <- incidence_i + new_incidence
    } else incidence_i <- incidence_i + colSums(new_incidence)
    # Compute the category of incidence
    incidence_0 <- incidence_1 <- incidence_2 <- incidence_i * 0
    incidence <- incidence_i / pop_mat[nrow(pop_mat),] * 1000000
    incidence_0[incidence <= 10] <- 1
    incidence_1[incidence > 10 & incidence <= thresh] <- 1
    incidence_2[incidence > thresh] <- 1
    # Compute the transmission potential
    if(bool_incidence == F | !is.null(reg_import)){
      FI_t <- (rbind(0 * data[date_obs >= time[i] - 50, rep(seq_len(94), n_sim)],
                     n_cases[date_ncas >= time[i] - 50 & date_ncas < time[i], ]) * 
                 mat_w[, rep(seq_len(94), n_sim)]) %>% colSums()
    } else{
      FI_t <- (rbind(data[date_obs >= time[i] - 50, rep(seq_len(94), n_sim)],
                     n_cases[date_ncas >= time[i] - 50 & date_ncas < time[i], ]) * 
                 mat_w[, rep(seq_len(94), n_sim)]) %>% colSums()      
    }
    
    # Compute the values of the predictors adding in the impact of incidence
    r0_ari <- mean_t$ar.exppred[i,] * 
      exp(incidence_1 * coef(model)["ar.cat1"]) * 
      exp(incidence_2 * coef(model)["ar.cat2"])
    r0_nei <- mean_t$ne.exppred[i,] * 
      exp(incidence_1 * coef(model)["ne.cat1"]) * 
      exp(incidence_2 * coef(model)["ne.cat2"])
    r0_eni <- mean_t$end.exppred[i,] * 
      exp(incidence_1 * coef(model)["end.cat1"]) * 
      exp(incidence_2 * coef(model)["end.cat2"])
    # Compute the average number of new cases per compartment
    new_ar_av <- r0_ari * FI_t
    new_en_av <- r0_eni 
    new_ne_av <- r0_nei * (sapply(seq_len(n_sim) - 1, function(X)
      (FI_t[seq_len(94) + X * 94] * wji) %>% colSums() %>% return) %>% c) 
    if(!is.null(reg_import)){
      new_en_av <- 0
    }
    # Generate the number of new cases per region
    n_cases[i,] <- rnbinom(length(new_ar_av),
                           mu = new_ar_av + new_ne_av + new_en_av, 
                           size = 1/(coef(model)["overdisp"]))
    if(!is.null(reg_import)){
      n_cases[i,][reg_import_nb] <- n_cases[i,][reg_import_nb] + 
        sum(dates_import == time[i])
    }
  }
  # Put each simulation set as an element of a list, and return the overall list
  sapply(seq_len(n_sim) - 1, function(X){
    return(list(n_cases[1:length(time),seq_len(94) + X * 94]))
  }) %>% return
  
}

# Function to integrate new conditions in the terms function"
term_predict <- function(model, conditions = list(), t = NULL, init_terms = NULL){
  # Set data$t to t
  if(!is.null(t)) model$control$data$t <- t
  # Create the initial terms output
  if(is.null(init_terms)) init_terms <- terms(model)
  init_terms$subset <- seq_len(init_terms$nTime)
  if(any(!is.element(names(conditions), init_terms$terms["name",])))
    stop("Names of conditions should match names of the model")
  # For each item in condition, set the value of the covariate in 
  # init_terms$terms to the new value
  for (i in names(conditions)){
    init_terms$terms["terms", init_terms$terms["name", ] == i] <-
      lapply(init_terms$terms["terms", init_terms$terms["name", ] == i],
             function(X){
               X[init_terms$subset, ] <- conditions[[i]]
               return(X)
             })
  }
  return(init_terms)
}

# Find coordinates of the regions of import
coordinates <- function(regions, types){
  # Import population centroid file, generated with function_centroid.R
  loc_centroid <- readRDS(file = "Data/location_centroids.RDS")
  # Keep only regions of interest
  loc_centroid <- loc_centroid[is.element(nuts3, regions),]
  setkey(loc_centroid, nuts3)
  loc_centroid <- loc_centroid[regions, reg := types]
  colnames(loc_centroid) <- c("long", "lat", "nuts3", "reg")
  setkey(loc_centroid, reg)
  return(loc_centroid)
}
