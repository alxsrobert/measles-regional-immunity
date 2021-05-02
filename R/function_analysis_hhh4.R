function_hhh4_daily <- function(pop_mat, area_mat, cov_mat, list_all_sim, 
                                distance_matrix, prop_gen1, mean_si, sd_si,
                                max_si, fun_wei, thresh = NA, start_coef = NULL){
  # Extract the number of simulations in list_all_sim
  n_sim <- length(list_all_sim$sim)
  # intialise list of all models
  list_models <- list()
  ## Run each model
  for(i in seq_len(n_sim)){
    data_i <- list_all_sim[[1]][[i]]
    # Generate the list of data needed to run the hhh4 models
    data_list <- prep_data(data = data_i, pop_mat = pop_mat, 
                           area_mat = area_mat, cov_mat = cov_mat, day = T, 
                           thresh = thresh, distance_matrix = distance_matrix, 
                           prop_gen1 = prop_gen1, mean_si = mean_si, 
                           sd_si = sd_si, max_si = max_si)
    ## Sts object 
    # observed is the potential for transmission, ie the convolution of the 
    # number of recent cases and the serial interval
    hhh4_sts <- sts(observed = data_list$previous, start = c(2009,1), 
                    frequency = 365, neighbourhood = data_list$distance_matrix)
    ## Control list
    hhh4_control <- list(
      end = list(f = (addSeason2formula(~1 + unvax + cat1 + cat2 + pop_mat
                                        , period = 365))),
      ar = list(f = addSeason2formula(~1 + unvax  + cat1 + cat2 + pop_mat + 
                                        area_mat, period = 365)),
      ne = list(f = addSeason2formula(~1 + unvax  + cat1 + cat2 + pop_mat 
                                      , period = 365),
                weights = fun_wei()), 
      data = list(unvax = log(data_list$unvax_ts), 
                  cat1 = data_list$cat_incidence1,
                  cat2 = data_list$cat_incidence2, 
                  pop = data_list$pop_mat,
                  pop_mat = log(data_list$pop_mat / 1000000), 
                  area_mat = log(data_list$area_mat),
                  response = data_i
      ), family = "NegBin1", verbose = T, start = list(fixed = start_coef))
    # Run the hhh4 model
    hhh4_run <- hhh4(stsObj = hhh4_sts, control = hhh4_control)
    # Add to the list of models
    list_models[[i]] <- hhh4_run
    
  }
  return(list_models)
}

function_hhh4_aggreg <- function(pop_mat, area_mat, cov_mat, list_all_sim, 
                                 distance_matrix, prop_gen1, fun_wei, 
                                 thresh = NA, len_agg = 10, start_coef = NULL){
  # Extract the number of simulations in list_all_sim
  n_sim <- length(list_all_sim$sim)
  # intialise list of all models
  list_models <- list()
  all_weeks <- as.Date(rownames(list_all_sim[[1]][[1]]))
  ## Run each model
  for(i in seq_len(n_sim)){
    data_i <- list_all_sim[[1]][[i]]
    ## Aggregate the data
    # Extract the first date of the dataset
    min_date <- min(rownames(data_i))
    # Vector of the dates of aggregation
    days_agg <- rep(all_weeks[all_weeks>= min_date], each = len_agg)
    days_agg <- days_agg[1:nrow(data_i)]
    # Initialise aggregated data set
    aggreg_data <- data_i
    rownames(aggreg_data) <- as.character(days_agg)
    # Aggregate the data
    aggreg_data <- t(sapply(by(aggreg_data,rownames(aggreg_data),colSums),
                            identity))
    
    # Generate the list of data needed to run the hhh4 models
    data_list <- prep_data(data = aggreg_data, area_mat = area_mat, 
                           pop_mat = pop_mat, cov_mat = cov_mat, day = F, 
                           thresh = thresh, distance_matrix = distance_matrix,
                           prop_gen1 = NULL, mean_si = NULL, 
                           sd_si = NULL, max_si = NULL)

    ## Sts object 
    # observed is the potential for transmission, ie the convolution of the 
    # number of recent cases and the serial interval
    hhh4_sts <- sts(observed = aggreg_data, start = c(2009,1), 
                    frequency = 365, neighbourhood = data_list$distance_matrix)
    ## Control list
    hhh4_control <- list(
      end = list(f = (addSeason2formula(~1 + unvax + cat1 + cat2 + pop_mat
                                        , period = 365))),
      ar = list(f = addSeason2formula(~1 + unvax  + cat1 + cat2 + pop_mat + 
                                        area_mat, period = 365)),
      ne = list(f = addSeason2formula(~1 + unvax  + cat1 + cat2 + pop_mat 
                                      , period = 365),
                weights = fun_wei()), 
      data = list(unvax = log(data_list$unvax_ts), 
                  cat1 = data_list$cat_incidence1,
                  cat2 = data_list$cat_incidence2, 
                  pop = data_list$pop_mat,
                  pop_mat = log(data_list$pop_mat / 1000000), 
                  area_mat = log(data_list$area_mat)
      ), family = "NegBin1", verbose = T, start = list(fixed = start_coef))
    
    # Run the hhh4 model
    hhh4_run <- hhh4(stsObj = hhh4_sts, control = hhh4_control)
    # Add to the list of models
    list_models[[i]] <- hhh4_run
    
  }
  return(list_models)
}