# Compute weights, dweights and d2weights. 
# Follows the structure of surveillance::W_powerlaw
W_exp_gravity_tot <- function(initial = c(1, 1)){
  weights <- function(params, ne_mat, data){
    # To speed up the computation, use the unique values of the population matrix
    uniq_pop_mat <- data$pop[!duplicated(data$pop),]
    # If there's only 1 unique value of the population matrix, 
    # re-format to matrix
    if(!is.matrix(uniq_pop_mat)){
      uniq_pop_mat <- matrix(uniq_pop_mat, nrow = 1)
      colnames(uniq_pop_mat) <- colnames(data$pop)
      rownames(uniq_pop_mat) <- rownames(data$pop)[1]
    }
    # extract the number of value of each unique item 
    dup <- cumsum(!duplicated(data$pop))
    # Extract the number of regions and the number of dates
    dim_uniq <- dim(uniq_pop_mat)[c(2,2,1)]
    nRegion <- dim_uniq[1]
    nTime <- dim_uniq[3]
    # Compute the population component of the gravity model
    pop_mat_params <- uniq_pop_mat^params[2] / uniq_pop_mat
    # Compute the distance component
    ne_mat_0 <- 1 * (ne_mat != 0)
    ne_mat_prod <- exp(-ne_mat*(params[1])) * ne_mat_0
    # Multiple the distance and the population components in a 3d matrix
    wji_uniq <- array( ne_mat_prod, dim = dim_uniq) *
      array(
        t(pop_mat_params)[rep(seq_len(nRegion), nRegion),],
        dim = dim_uniq)
    # Repeat by the number of duplicates
    wji_uniq[,,dup]
  }
  dweights <- function(params, ne_mat, data){
    # Compute first order derivatives 
    uniq_pop_mat <- data$pop[!duplicated(data$pop),]
    if(!is.matrix(uniq_pop_mat)){
      uniq_pop_mat <- matrix(uniq_pop_mat, nrow = 1)
      colnames(uniq_pop_mat) <- colnames(data$pop)
      rownames(uniq_pop_mat) <- rownames(data$pop)[1]
    }
    dup <- cumsum(!duplicated(data$pop))
    dim_uniq <- dim(uniq_pop_mat)[c(2,2,1)]
    nRegion <- dim_uniq[1]
    nTime <- dim_uniq[3]
    
    pop_mat_params <- uniq_pop_mat^params[2] / uniq_pop_mat
    ne_mat_0 <- 1 * (ne_mat != 0)
    
    
    ne_mat_prod <- exp(-ne_mat*(params[1])) * ne_mat_0
    
    wji_uniq <- array( ne_mat_prod, dim = dim_uniq) *
      array(
        t(pop_mat_params)[rep(seq_len(nRegion), nRegion),],
        dim = dim_uniq)
    
    logpop <- array(t(log(uniq_pop_mat))[rep(seq_len(nRegion), nRegion),],
                    dim = dim_uniq)
    array_ne_mat <- array(ne_mat, dim = dim_uniq)
    
    deriv_1 <- - array_ne_mat * wji_uniq
    
    deriv_2 <- logpop * wji_uniq
    
    deriv_1[is.na(deriv_1)] <- 0
    deriv_2[is.na(deriv_2)] <- 0
    list(deriv_1[,,dup], deriv_2[,,dup])
  }
  d2weights <- function(params, ne_mat, data){
    # Compute second order derivatives 
    uniq_pop_mat <- data$pop[!duplicated(data$pop),]
    if(!is.matrix(uniq_pop_mat)){
      uniq_pop_mat <- matrix(uniq_pop_mat, nrow = 1)
      colnames(uniq_pop_mat) <- colnames(data$pop)
      rownames(uniq_pop_mat) <- rownames(data$pop)[1]
    }
    dup <- cumsum(!duplicated(data$pop))
    dim_uniq <- dim(uniq_pop_mat)[c(2,2,1)]
    nRegion <- dim_uniq[1]
    nTime <- dim_uniq[3]
    
    pop_mat_params <- uniq_pop_mat^params[2] / uniq_pop_mat
    ne_mat_0 <- 1 * (ne_mat != 0)
    
    ne_mat_prod <- exp(-ne_mat*(params[1])) * ne_mat_0 
    
    wji_uniq <- array( ne_mat_prod, dim = dim_uniq) *
      array(
        t(pop_mat_params)[rep(seq_len(nRegion), nRegion),],
        dim = dim_uniq)
    
    logpop <- array(t(log(uniq_pop_mat))[rep(seq_len(nRegion), nRegion),],
                    dim = dim_uniq)
    array_ne_mat <- array(ne_mat, dim = dim_uniq)
    
    
    dderiv_1 <- - array_ne_mat * - array_ne_mat * wji_uniq
    dderiv_2 <- logpop * - array_ne_mat * wji_uniq
    dderiv_3 <- logpop * logpop * wji_uniq
    
    dderiv_1[is.na(dderiv_1)] <- 0
    dderiv_2[is.na(dderiv_2)] <- 0
    dderiv_3[is.na(dderiv_3)] <- 0
    list(dderiv_1[,,dup], dderiv_2[,,dup], dderiv_3[,,dup])
    
  }
  list(w=weights, dw=dweights, d2w=d2weights, initial=initial)
}
# No need for distance parameter if only neighbours are connected
W_exp_gravity_nei <- function(initial = c(1)){
  weights <- function(params, ne_mat, data){
    # To speed up the computation, use the unique values of the population matrix
    uniq_pop_mat <- data$pop[!duplicated(data$pop),]
    # If there's only 1 unique value of the population matrix, 
    # re-format to matrix
    if(!is.matrix(uniq_pop_mat)){
      uniq_pop_mat <- matrix(uniq_pop_mat, nrow = 1)
      colnames(uniq_pop_mat) <- colnames(data$pop)
      rownames(uniq_pop_mat) <- rownames(data$pop)[1]
    }
    # extract the number of value of each unique item 
    dup <- cumsum(!duplicated(data$pop))
    # Extract the number of regions and the number of dates
    dim_uniq <- dim(uniq_pop_mat)[c(2,2,1)]
    nRegion <- dim_uniq[1]
    nTime <- dim_uniq[3]
    # Compute the population component of the gravity model
    pop_mat_params <- uniq_pop_mat^params[1] / uniq_pop_mat
    # Compute the distance component
    ne_mat_0 <- 1 * (ne_mat != 0)
    ne_mat_prod <- exp(-ne_mat*(1)) * ne_mat_0
    # Multiple the distance and the population components in a 3d matrix
    wji_uniq <- array( ne_mat_prod, dim = dim_uniq) *
      array(
        t(pop_mat_params)[rep(seq_len(nRegion), nRegion),],
        dim = dim_uniq)
    # Repeat by the number of duplicates
    wji_uniq[,,dup]
  }
  dweights <- function(params, ne_mat, data){
    # Compute first order derivatives 
    uniq_pop_mat <- data$pop[!duplicated(data$pop),]
    if(!is.matrix(uniq_pop_mat)){
      uniq_pop_mat <- matrix(uniq_pop_mat, nrow = 1)
      colnames(uniq_pop_mat) <- colnames(data$pop)
      rownames(uniq_pop_mat) <- rownames(data$pop)[1]
    }
    dup <- cumsum(!duplicated(data$pop))
    dim_uniq <- dim(uniq_pop_mat)[c(2,2,1)]
    nRegion <- dim_uniq[1]
    nTime <- dim_uniq[3]
    
    pop_mat_params <- uniq_pop_mat^params[1] / uniq_pop_mat
    ne_mat_0 <- 1 * (ne_mat != 0)
    
    
    ne_mat_prod <- exp(-ne_mat*(1)) * ne_mat_0
    
    wji_uniq <- array( ne_mat_prod, dim = dim_uniq) *
      array(
        t(pop_mat_params)[rep(seq_len(nRegion), nRegion),],
        dim = dim_uniq)
    
    logpop <- array(t(log(uniq_pop_mat))[rep(seq_len(nRegion), nRegion),],
                    dim = dim_uniq)
    array_ne_mat <- array(ne_mat, dim = dim_uniq)
    
    deriv_2 <- logpop * wji_uniq
    
    deriv_2[is.na(deriv_2)] <- 0
    list(deriv_2[,,dup])
  }
  d2weights <- function(params, ne_mat, data){
    # Compute second order derivatives 
    uniq_pop_mat <- data$pop[!duplicated(data$pop),]
    if(!is.matrix(uniq_pop_mat)){
      uniq_pop_mat <- matrix(uniq_pop_mat, nrow = 1)
      colnames(uniq_pop_mat) <- colnames(data$pop)
      rownames(uniq_pop_mat) <- rownames(data$pop)[1]
    }
    dup <- cumsum(!duplicated(data$pop))
    dim_uniq <- dim(uniq_pop_mat)[c(2,2,1)]
    nRegion <- dim_uniq[1]
    nTime <- dim_uniq[3]
    
    pop_mat_params <- uniq_pop_mat^params[1] / uniq_pop_mat
    ne_mat_0 <- 1 * (ne_mat != 0)
    
    ne_mat_prod <- exp(-ne_mat*(1)) * ne_mat_0 
    
    wji_uniq <- array( ne_mat_prod, dim = dim_uniq) *
      array(
        t(pop_mat_params)[rep(seq_len(nRegion), nRegion),],
        dim = dim_uniq)
    
    logpop <- array(t(log(uniq_pop_mat))[rep(seq_len(nRegion), nRegion),],
                    dim = dim_uniq)
    array_ne_mat <- array(ne_mat, dim = dim_uniq)
    
    
    dderiv_3 <- logpop * logpop * wji_uniq
    
    dderiv_3[is.na(dderiv_3)] <- 0
    list(dderiv_3[,,dup])
    
  }
  list(w=weights, dw=dweights, d2w=d2weights, initial=initial)
}