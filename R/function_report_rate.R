report_rate <- function(report_rate, data){
  # Initialise matrix of reported cases
  data_reported <- data * 0
  # Compute the number of reported cases 
  n_reported_cases <- round(sum(data) * report_rate)
  # Vector of the coordinate of each case in the data
  all_cases <- rep(seq_along(data), as.numeric(data))
  # Select a proportion of report_rate cases, and add them to data_reported  
  report_cases <- sort(sample(x = all_cases, n_reported_cases))
  data_reported[unique(report_cases)] <- as.numeric(table(report_cases)) 
  
  return(data_reported)
}
