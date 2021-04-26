prep_data <- function(data, pop_mat, area_mat, cov_mat, day, thresh, 
                      distance_matrix, prop_gen1, mean_si, sd_si, max_si){
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
  
  if(is.na(thresh)) thresh <- quantile(incidence, 2/3)
  print(thresh)
  # Compute the three categories of incidence
  cat_incidence1 <- cat_incidence2 <- incidence * 0
  cat_incidence1[incidence > 10 & incidence <= thresh] <- 1
  cat_incidence2[incidence >= thresh] <- 1
  
  # Order pop area and coverage matrices
  pop_mat_hhh4 <- pop_mat[rownames(data), colnames(data)]
  area_mat_hhh4 <- area_mat[rownames(data), colnames(data)]
  distance_matrix_hhh4 <- distance_matrix[colnames(data), colnames(data)]
  cov_mat_hhh4 <- (100 - cov_mat[rownames(data), colnames(data)])/100
  
  ## Compute the potential for infection (convolve SI with nb cases)
  if(day == TRUE){
    # Compute the serial interval
    w_dens <- serial_interval(prop_gen1 = prop_gen1, mean_w1 = mean_si, 
                              sd_w1 = sd_si, max_days = max_si)
    # initialise n_Cases, the matrix containing the potential for infections
    n_cases <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    t_min <- min(as.Date(rownames(data)))
    
    ## Convolve the number of recent cases per region with the serial interval
    for (t in seq_len(nrow(data))[-1]){
      if(t<=50) t_i_min <- t_min else 
        t_i_min <- as.Date(rownames(data)[t - 50])
      date_t <- as.Date(rownames(data)[t])
      dates_t <- as.Date(t_i_min:(date_t-1),  origin = "1970-01-01")
      cases_tot <- data[as.character(dates_t),] * 
        matrix(w_dens[length(dates_t):1], 
               nrow = length(dates_t), ncol = ncol(data), byrow = F)
      
      n_cases[t-1,] <- colSums(cases_tot)
      colnames(n_cases) <- colnames(data)
    }
  } else n_cases <- NULL
    
  
  return(list(pop_mat = pop_mat_hhh4, area_mat = area_mat_hhh4,
              cat_incidence1 = cat_incidence1, cat_incidence2 = cat_incidence2, 
              distance_matrix = distance_matrix_hhh4, 
              unvax_ts = cov_mat_hhh4, previous = n_cases))
}
