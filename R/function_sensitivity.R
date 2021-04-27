sens_SI <- function(vals, list_all_sim, pop_mat, area_mat, cov_mat, fun_wei,
                    distance_matrix, mean_si, sd_si, max_si, thresh = NA,
                    start = NULL){
  # Run the model, with the proportion of transmission without missing cases
  # indicated in the vector vals
  for(i in seq_along(vals)){
    model_i <- function_hhh4_daily(list_all_sim = list_all_sim, 
                                   pop_mat = pop_mat, area_mat = area_mat, 
                                   cov_mat = cov_mat, thresh = thresh, 
                                   distance_matrix = distance_matrix,
                                   prop_gen1 = vals[i], mean_si = mean_si, 
                                   sd_si = sd_si, max_si = max_si,
                                   fun_wei = fun_wei, start = start)[[1]]
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
                     start = NULL, prop_gen1){
  
  for(i in seq_len(n)){
    cov_mat <- importation_coverage(corres = corres, len_break = 1, rd_cov = T)
    model_i <- function_hhh4_daily(list_all_sim = list_all_sim, 
                                   pop_mat = pop_mat, area_mat = area_mat, 
                                   cov_mat = cov_mat, thresh = thresh, 
                                   distance_matrix = distance_matrix,
                                   prop_gen1 = prop_gen1, mean_si = mean_si, 
                                   sd_si = sd_si, max_si = max_si,
                                   fun_wei = fun_wei, start = start)[[1]]
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
