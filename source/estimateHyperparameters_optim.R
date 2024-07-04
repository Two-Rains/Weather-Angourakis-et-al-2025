source("source/weatherModel.R")

estimateHyperparameters_optim <- function(curves, objectiveFunction, method, lower, upper, initialGuess)
{
  bestEstimationCurves <- list()
  bestEstimationFits <- list()
  
  for (i in 1:length(curves))
  {
    observedData <- curves[[i]]
    
    # Use the least squares method to estimate the parameter values
    
    fit <- optim(initialGuess, 
                 objectiveFunction,
                 observedData = observedData,
                 method = method, 
                 lower = lower, 
                 upper = upper)
    
    bestEstimationFits[[i]] <- fit
  }
  
  estimated_summary <- calculate_mean_sd(get_list_par_from_fit(bestEstimationFits))
  
  return(list(curve_fits = bestEstimationFits, summary = estimated_summary))
}

calculate_mean_sd <- function(list_of_vectors) {
  # Initialize empty vectors to store mean and standard deviation
  mean_vector <- numeric(length(list_of_vectors[[1]]))
  sd_vector <- numeric(length(list_of_vectors[[1]]))
  
  # Calculate mean and standard deviation for each position
  for (i in 1:length(list_of_vectors[[1]])) {
    values <- sapply(list_of_vectors, function(x) x[i])
    mean_vector[i] <- mean(values)  # Calculate mean
    sd_vector[i] <- sd(values)  # Calculate standard deviation
  }
  
  return(list(mean = mean_vector, sd = sd_vector))  # Return a list of mean and sd vectors
}

get_list_par_from_fit <- function(listOffitObjects) {
  
  par_list <- list()
  for (i in 1:length(listOffitObjects))
  {
    par_list[[i]] <- listOffitObjects[[i]]$par
  }
  return(par_list)
}
