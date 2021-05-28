# Generate the calibration scores for a list of models
calibration_model <- function(model_list, daily, all_dates,
                              w_dens, period_pred = 10, k = 1000,
                              new_fit = FALSE){
  # Compute the number of model
  n_model <- length(model_list)
  # Compute the number of regions
  n_reg <- model_list[[1]]$nUnit
  # Initialise an array containing all samples for each day of the study period
  k_sample <- array(dim = c(period_pred, n_reg, k))
  # Matrix replicating the serial interval for each region 
  mat_w <- matrix(rev(w_dens), nrow = length(w_dens), ncol = n_reg, byrow = F)
  # Initialise bias, sharpness, rps and logs matrices
  mat_bias <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_shar <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_rps <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_log <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_tot <- array(dim = c(length(all_dates), 5, n_model))
  rownames(mat_tot) <- all_dates
  # Array with the proportion of simulations where the number of cases is 
  # equal to the number of cases in the observation
  px <- array(dim = c(length(all_dates), n_reg, n_model))
  # Array with the proportion of simulations where the number of cases is 
  # equal to the number of cases in the observation - 1
  pxm1 <- array(dim = c(length(all_dates), n_reg, n_model))
  # Calibration for each model in model_list
  for (i in seq_len(n_model)){
    model <- model_list[[i]]
    data <- model$control$data$response
    # If daily model, run simulations
    if(daily){
      fixed_coef <- fixef(model)
      for(j in seq_along(all_dates)){
        # Extract the date at i
        date <- all_dates[j]
        # Extract the input parameters of the hhh4 model
        control <- model$control
        stsobj <- model$stsObj
        obs <- stsobj@observed
        # Set the last observation to "date"
        max_date_obs <- rownames(model$data$response)[date]
        data_j <- data[date + 1:period_pred,] %>% colSums()
        control$subset <- 2:date
        new_obs <- obs[1:max(control$subset),]
        new_resp <- control$data$response[1:max(control$subset),]
        # if new_fit = T, fit the model at date 
        if(new_fit){
          # Use the coefficients of the last model as starting points to 
          # speed up the fitting procedure 
          control$start <- list(fixed = fixed_coef)
          control$verbose <- FALSE
          model_j <- hhh4(stsObj = stsobj, control = control)
          fixed_coef <- fixef(model_j)
        } else model_j <- model
        # Simulate k samples over period_pred
        k_sum <- simulate_calib(model_j = model_j, period_pred = period_pred, 
                                observed = obs, k_sample = k_sample, 
                                control = control, j = date, mat_w = mat_w)
        # If any of the data point was not reached, run another k samples
        while(any(rowSums(k_sum == data_j) == 0)){
          k_sum2 <- simulate_calib(model_j = model_j, period_pred = period_pred, 
                                   observed = obs, k_sample = k_sample, 
                                   control = control, j = date, mat_w = mat_w)
          k_sum <- cbind(k_sum, k_sum2)
        }
        mat_tot[j, ,i] <- quantile(colSums(k_sum), 
                                   probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
        # Compute the logs score for each entry
        mat_log[j,,i] <- -log(rowSums(k_sum == data_j)/ncol(k_sum))
        # For each region, compute the bias, sharpness, and RPS score (scoringutils)
        for(reg in seq_len(n_reg)){
          mat_pred <- matrix(k_sum[reg,], nrow = 1)
          mat_bias[j,reg, i] <- bias(data_j[reg], mat_pred)
          mat_shar[j,reg, i] <- sharpness(mat_pred)
          mat_rps[j,reg, i] <- crps(data_j[reg], mat_pred)
          # Compute px and pxm1 to generate the PIT histograms
          px[j, reg, i] <- sum(mat_pred <= data_j[reg])/length(mat_pred)
          if(data_j[reg] == 0) pxm1[j, reg, i] <- 0 else
            pxm1[j, reg, i] <- sum(mat_pred <= data_j[reg]-1)/length(mat_pred)
        }
      }
    } else {
      # If aggregated model, then use the built-in function
      list_pred <- oneStepAhead(model, tp = c(all_dates[1], max(all_dates)), 
                                type = "rolling", which.start = "final")
      pred <- list_pred$pred
      psi <- list_pred$psi
      obs <- list_pred$obs
      # For each region, compute the bias, sharpness, and RPS score (scoringutils)
      for(j in seq_along(all_dates)){
        for(reg in seq_len(n_reg)){
          mat_pred <- matrix(rnbinom(n = k, mu = pred[j, reg], 
                                     size = exp(psi[j,])), nrow = 1)
          
          mat_bias[j,reg, i] <- bias(obs[j, reg], mat_pred)
          mat_shar[j,reg, i] <- sharpness(mat_pred)
          mat_rps[j,reg, i] <- crps(obs[j, reg], mat_pred)
          mat_log[j,,i] <- -log(((mat_pred == obs[j, reg]) %>% rowSums())/k)
          px[j, reg, i] <- sum(mat_pred <= obs[j,reg])/length(mat_pred)
          
          if(obs[j,reg] == 0) pxm1[j, reg, i] <- 0 else
            pxm1[j, reg, i] <- sum(mat_pred <= obs[j,reg]-1)/length(mat_pred)
        }
      }
      
    }
  }
  return(list(bias = mat_bias, sharpness = mat_shar, rps = mat_rps, 
              log = mat_log, px = px, pxm1 = pxm1, tot = mat_tot))
}

# Simulate over the calibration period
simulate_calib <- function(model_j, period_pred, observed, k_sample, control, 
                           j, mat_w){
  # k: number of samples
  k <- dim(k_sample)[3]
  ## Compute the regional coefficients at this time step
  term_j <- terms(model_j)
  term_j$subset <- j + 1:period_pred
  hhh_j <- meanHHH(coef(model_j), term_j)
  # Compute the connectivity matrix
  weight_j <- getNEweights(model_j)[,,j]
  # Transmission potential at j
  last_obs <-  observed[j,]
  # Overdispersion parameter
  overdisp <- 1/coef(model_j)["overdisp"]
  
  for (t in seq_len(period_pred)){
    if(t == 1){
      # Compute the average number of new cases
      n_cases <- (hhh_j$ar.exppred[1,] * last_obs + 
                    hhh_j$end.exppred[1,] + 
                    hhh_j$ne.exppred[1,] * ((last_obs * weight_j) %>% 
                                              colSums()))
      # Draw the number of new cases in each sample
      k_sample[t,,] <- t(sapply(n_cases, function(X)
        rnbinom(k, mu = X, size = overdisp)))
    }else{
      ## Compute the transmission potential
      # Compute the part that comes from the data
      data_nb_fit <- control$data$response[
        j - (length(w_dens):1 - 1),][-c(1:(t-1)),]
      pot_transmission <- colSums(data_nb_fit * mat_w[-(nrow(mat_w) - 0:(t-2)),])
      # Compute the part specific to each sample (i.e. transmission since t ==1)
      for(iter in seq_len(t-1)){
        pot_transmission <- pot_transmission + 
          k_sample[iter,,] * w_dens[t - iter]
      }
      # Compute the number of samples with only 0s
      if(t == 2) {
        k_0 <- which(((k_sample[(1:t-1),,]) %>% colSums()) == 0)
      } else{
        k_0 <- which(((k_sample[(1:t-1),,]) %>% colSums() %>% colSums()) == 0)
      }
      ## Compute transmission potential * weights:
      # If more than 10% of the samples are only 0s, then only calculate for 
      # the samples with at least 1 case (to accelerate this section)
      if(length(k_0)>(k / 10)){
        uniq_pot <- pot_transmission[, -k_0[-1]]
        uniq_weight <- apply(uniq_pot, 2, function(X)
          return(colSums(X * weight_j)))
        weight_t <- pot_transmission * 0
        weight_t[,k_0] <- uniq_weight[,k_0[1]]
        weight_t[,-k_0] <- uniq_weight[,-k_0[1]]
      } else{
        # Otherwise: multiple pot_transmission in each region to the 
        # connectivity matrix
        weight_t <- apply(pot_transmission, 2, function(X) 
          return(colSums(X * weight_j)))
      }
      # Compute the average number of new cases
      n_cases <- (hhh_j$ar.exppred[t,] * pot_transmission + 
                    hhh_j$end.exppred[t,] + 
                    hhh_j$ne.exppred[t,] * weight_t)
      vec_cases <- c(n_cases)
      # Draw the number of new cases in each sample
      k_sample[t,,] <- t(matrix(rnbinom(vec_cases,
                                        mu = vec_cases,size = overdisp),
                                nrow = k, byrow = T))
    }
  }
  # compute the number of cases over period_pred in each sample
  k_sum <- colSums(k_sample)
  return(k_sum)
}

