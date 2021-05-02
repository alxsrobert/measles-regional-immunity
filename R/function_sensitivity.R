sens_SI <- function(vals, list_all_sim, pop_mat, area_mat, cov_mat, fun_wei,
                    distance_matrix, mean_si, sd_si, max_si, thresh = NA,
                    start_coef = NULL){
  # Run the model, with the proportion of transmission without missing cases
  # indicated in the vector vals
  for(i in seq_along(vals)){
    model_i <- function_hhh4_daily(list_all_sim = list_all_sim, 
                                   pop_mat = pop_mat, area_mat = area_mat, 
                                   cov_mat = cov_mat, thresh = thresh, 
                                   distance_matrix = distance_matrix,
                                   prop_gen1 = vals[i], mean_si = mean_si, 
                                   sd_si = sd_si, max_si = max_si,
                                   fun_wei = fun_wei, start_coef = start_coef)[[1]]
    # At the first iteration, initialise the mean coefficient matrix, 
    # and the 95% CI matrices
    if(i == 1){
      coef_mod <- matrix(ncol = length(coef(model_i)), nrow = length(vals))
      max_coef <- matrix(ncol = length(coef(model_i)), nrow = length(vals))
      min_coef <- matrix(ncol = length(coef(model_i)), nrow = length(vals))
      colnames(coef_mod) <- colnames(max_coef) <- colnames(min_coef) <- 
        names(coef(model_i))
    }
    coef_mod[i, ] <- coef(model_i)
    min_coef[i,] <- confint(model_i)[,1]
    max_coef[i,] <- confint(model_i)[,2]
  }
  return(list(coefs = coef_mod, min_coef = min_coef, max_coef = max_coef))
}

sens_vax <- function(n, list_all_sim, pop_mat, area_mat, fun_wei, corres,
                     distance_matrix, mean_si, sd_si, max_si, thresh = NA,
                     start_coef = NULL, prop_gen1){
  
  for(i in seq_len(n)){
    cov_mat <- importation_coverage(corres = corres, len_break = 1, rd_cov = T)
    model_i <- function_hhh4_daily(list_all_sim = list_all_sim, 
                                   pop_mat = pop_mat, area_mat = area_mat, 
                                   cov_mat = cov_mat, thresh = thresh, 
                                   distance_matrix = distance_matrix,
                                   prop_gen1 = prop_gen1, mean_si = mean_si, 
                                   sd_si = sd_si, max_si = max_si,
                                   fun_wei = fun_wei, start_coef = start_coef)[[1]]
    if(i == 1){
      coef_mod <- matrix(ncol = length(coef(model_i)), nrow = n)
      max_coef <- matrix(ncol = length(coef(model_i)), nrow = n)
      min_coef <- matrix(ncol = length(coef(model_i)), nrow = n)
      colnames(coef_mod) <- colnames(max_coef) <- colnames(min_coef) <- 
        names(coef(model_i))
    }
    coef_mod[i, ] <- coef(model_i)
    min_coef[i,] <- confint(model_i)[,1]
    max_coef[i,] <- confint(model_i)[,2]
  }
  return(list(coefs = coef_mod, min_coef = min_coef, max_coef = max_coef))
}

sens_weekday <- function(list_all_sim, pop_mat, area_mat, fun_wei, prop_gen1,
                         distance_matrix, mean_si, sd_si, max_si, cov_mat,
                         thresh = NA, start_coef = NULL){
  data <- list_all_sim[[1]][[1]]
  # Generate the list of data needed to run the hhh4 models
  data_list <- prep_data(data = data, pop_mat = pop_mat, 
                         area_mat = area_mat, cov_mat = cov_mat, day = T, 
                         thresh = thresh, distance_matrix = distance_matrix, 
                         prop_gen1 = prop_gen1, mean_si = mean_si, 
                         sd_si = sd_si, max_si = max_si)
  ## Compute the day of the week matrix
  weekday_ts <- data * 0
  weekday_ts[is.element(weekdays(as.Date(rownames(weekday_ts))), 
                        c("Saturday", "Sunday"))] <- 1
  
  ## Sts object 
  # observed is the potential for transmission, ie the convolution of the 
  # number of recent cases and the serial interval
  hhh4_sts <- sts(observed = data_list$previous, start = c(2009,1), 
                  frequency = 365, neighbourhood = data_list$distance_matrix)
  ## Control list
  hhh4_control <- list(
    end = list(f = (addSeason2formula(~1 + unvax + cat1 + cat2 + pop_mat + weekday
                                      , period = 365))),
    ar = list(f = addSeason2formula(~1 + unvax  + cat1 + cat2 + pop_mat + 
                                      area_mat + weekday, period = 365)),
    ne = list(f = addSeason2formula(~1 + unvax  + cat1 + cat2 + pop_mat  + weekday 
                                    , period = 365),
              weights = fun_wei()), 
    data = list(unvax = log(data_list$unvax_ts), 
                cat1 = data_list$cat_incidence1,
                cat2 = data_list$cat_incidence2, 
                weekday = weekday_ts,
                pop = data_list$pop_mat,
                pop_mat = log(data_list$pop_mat / 1000000), 
                area_mat = log(data_list$area_mat),
                response = data
    ), family = "NegBin1", verbose = T, start = list(fixed = start_coef))
  # Run the hhh4 model
  hhh4_run <- hhh4(stsObj = hhh4_sts, control = hhh4_control)
  return(hhh4_run)
}