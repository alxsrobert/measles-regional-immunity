conv <- function(w_dens, max_days){
  # Compute log serial interval
  x <- log(w_dens)
  y <- rev(log(w_dens))
  # Compute the length of the serial interval
  n <- m <- length(x)
  
  r <- lapply(1:(n+m-1), function(k){
    i <- seq_len(n)
    i <- i[(k > (n - i))]
    w <- x[k-m+i]+y[i]
    
    total <- w[1]
    
    if (length(w)<2) return(total)
    
    for (i in 2:length(w)){
      total <- max(total, w[i]) + log(exp(total - max(total, w[i])) + 
                                        exp(w[i] - max(total, w[i])))
    }
    return(total)
  })[1 + 1:max_days] %>% unlist %>% exp
  return(r)
}

serial_interval <- function(prop_gen1, mean_w1, sd_w1, max_days = 50){
  # Generate the distribution of the serial interval (no missing case)
  w1 <- dnorm(x = 0:max_days, mean = mean_w1, sd = sd_w1)
  # Generate the distribution of the serial interval (missing ancestor)
  w0 <- dnorm(x = 1:max_days, mean = 0, sd = sqrt(2 * sd_w1**2))
  # Normalise distribution (>= 1 day)
  w0 <- w0/sum(w0)
  # Generate the distribution of the serial interval (missing generation)
  w2 <- conv(w1, max_days = max_days)
  w1 <- w1[-1]
  # Compute the composite serial interval
  w_dens <- w0 * ((1-prop_gen1)/2) + w1 * prop_gen1 + w2 * ((1-prop_gen1)/2)
  return(w_dens)
}