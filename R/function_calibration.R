calibration_model <- function(model_list, daily, all_dates,
                              w_dens, period_pred = 10, k = 1000,
                              new_fit = FALSE){

  n_reg <- model_list[[1]]$nUnit
  n_model <- length(model_list)
  
  k_sample <- array(dim = c(period_pred, n_reg, k))
  mat_w <- matrix(rev(w_dens), nrow = length(w_dens), ncol = n_reg, byrow = F)
  
  mat_bias <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_shar <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_rps <- array(dim = c(length(all_dates), n_reg, n_model))
  mat_log <- array(dim = c(length(all_dates), n_reg, n_model))
  px <- array(dim = c(length(all_dates), n_reg, n_model))
  pxm1 <- array(dim = c(length(all_dates), n_reg, n_model))
  
  for (i in seq_len(n_model)){
    model <- model_list[[i]]
    data <- model$control$data$response
    if(daily){
      for(j in seq_along(all_dates)){
        print(j)
        date <- all_dates[j]
        control <- model$control
        stsobj <- model$stsObj
        
        obs <- stsobj@observed
        
        max_date_obs <- rownames(model$data$response)[date]
        
        data_j <- data[date + 1:period_pred,] %>% colSums()
        control$subset <- 2:date
        new_obs <- obs[1:max(control$subset),]
        new_resp <- control$data$response[
          1:max(control$subset),]
        if(new_fit){
          control$start <- list(fixed = fixef(model))
          model_j <- hhh4(stsObj = stsobj, control = control)
        } else model_j <- model
        
        k_sum <- simulate_calib(model_j = model_j, period_pred = period_pred, 
                                observed = obs, k_sample = k_sample, 
                                control = control, j = date, mat_w = mat_w)
        while(any(rowSums(k_sum == data_j) == 0)){
          k_sum2 <- simulate_calib(model_j = model_j, period_pred = period_pred, 
                                   observed = obs, k_sample = k_sample, 
                                   control = control, j = date, mat_w = mat_w)
          k_sum <- cbind(k_sum, k_sum2)
        }
        mat_log[j,,i] <- -log(rowSums(k_sum == data_j)/ncol(k_sum))
        for(reg in seq_len(n_reg)){
          mat_pred <- matrix(k_sum[reg,], nrow = 1)
          mat_bias[j,reg, i] <- bias(data_j[reg], mat_pred)
          mat_shar[j,reg, i] <- sharpness(mat_pred)
          mat_rps[j,reg, i] <- crps(data_j[reg], mat_pred)
          
          px[j, reg, i] <- sum(mat_pred <= data_j[reg])/length(mat_pred)
          if(data_j[reg] == 0) pxm1[j, reg, i] <- 0 else
            pxm1[j, reg, i] <- sum(mat_pred <= data_j[reg]-1)/length(mat_pred)
        }
      }
    } else {
      list_pred <- oneStepAhead(model, tp = c(all_dates[1], max(all_dates)), 
                                type = "rolling", which.start = "final")
      pred <- list_pred$pred
      psi <- list_pred$psi
      obs <- list_pred$obs
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
              log = mat_log, px = px, pxm1 = pxm1))
}


simulate_calib <- function(model_j, period_pred, observed, k_sample, control, 
                           j, mat_w){
  k <- dim(k_sample)[3]
  
  term_j <- terms(model_j)
  term_j$subset <- j + 1:period_pred
  
  weight_j <- getNEweights(model_j)[,,j]
  
  hhh_j <- meanHHH(coef(model_j), term_j)
  
  last_obs <-  observed[j,]
  
  overdisp <- 1/coef(model_j)["overdisp"]
  
  for (t in seq_len(period_pred)){
    if(t == 1){
      n_cases <- (hhh_j$ar.exppred[1,] * last_obs + 
                    hhh_j$end.exppred[1,] + 
                    hhh_j$ne.exppred[1,] * ((last_obs * weight_j) %>% 
                                              colSums()))
      k_sample[t,,] <- t(sapply(n_cases, function(X)
        rnbinom(k, mu = X, size = overdisp)))
    }else{
      data_nb_fit <- control$data$response[
        j - (length(w_dens):1 - 1),][-c(1:(t-1)),]
      pot_transmission <- colSums(data_nb_fit * mat_w[-(nrow(mat_w) - 0:(t-2)),])
      for(iter in seq_len(t-1)){
        pot_transmission <- pot_transmission + 
          k_sample[iter,,] * w_dens[t - iter]
      }
      if(t == 2) {
        k_0 <- which(((k_sample[(1:t-1),,]) %>% colSums()) == 0)
      } else{
        k_0 <- which(((k_sample[(1:t-1),,]) %>% colSums() %>% colSums()) == 0)
      }
      if(length(k_0)>(k / 10)){
        uniq_pot <- pot_transmission[, -k_0[-1]]
        uniq_weight <- apply(uniq_pot, 2, function(X)
          return(colSums(X * weight_j)))
        weight_t <- pot_transmission * 0
        weight_t[,k_0] <- uniq_weight[,k_0[1]]
        weight_t[,-k_0] <- uniq_weight[,-k_0[1]]
      } else{
        weight_t <- apply(pot_transmission, 2, function(X) 
          return(colSums(X * weight_j)))
      }
      n_cases <- (hhh_j$ar.exppred[t,] * pot_transmission + 
                    hhh_j$end.exppred[t,] + 
                    hhh_j$ne.exppred[t,] * weight_t)
      vec_cases <- c(n_cases)
      k_sample[t,,] <- t(matrix(rnbinom(vec_cases,
                                        mu = vec_cases,size = overdisp),
                                nrow = k, byrow = T))
    }
  }
  k_sum <- colSums(k_sample)
  return(k_sum)
}

