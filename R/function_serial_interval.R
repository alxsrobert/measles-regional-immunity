conv <- function(w_dens, max_days){
  x <- log(w_dens)
  y <- rev(log(w_dens))
  
  n <- length(x)
  m <- length(y)
  
  r <- lapply(1:(n+m-1), function(k){
    i <- 1:max(m,n)
    i <- i[((i<=m) & ((k-m+i) <= n)) & ((k-m+i) > 0)]
    w <- x[k-m+i]+y[i]
    
    total <- w[1]
    
    if (length(w)<2) return(total)
    
    for (i in 2:length(w)){
      total <- max(total, w[i]) + log(exp(total - max(total, w[i])) + 
                                        exp(w[i] - max(total, w[i])))
    }
    return(total)
  })[1:max_days] %>% unlist %>% exp
  return(r[seq_along(w_dens)])
}

serial_interval <- function(prop_gen1, mean_w1, sd_w1, max_days = 50){
  w1 <- dnorm(x = 1:max_days, mean = mean_w1, sd = sd_w1)
  w0 <- dnorm(x = 1:max_days, mean = 0, sd = sqrt(2 * sd_w1**2))
  w0 <- w0/sum(w0)
  
  w2 <- conv(w1, max_days = max_days)
  
  w_dens <- w0 * ((1-prop_gen1)/2) + w1 * prop_gen1 + w2 * ((1-prop_gen1)/2)
  return(w_dens)
}