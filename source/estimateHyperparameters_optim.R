source("source/weatherModel.R")

estimate_hyperparameters_optim <- function(curves, objective_function, method, lower, upper, initial_guess) {
  best_estimation_fits <- lapply(curves, function(observed_data) {
    optim(initial_guess, 
          objective_function,
          observedData = observed_data,
          method = method, 
          lower = lower, 
          upper = upper)
  })
  
  estimated_summary <- calculate_mean_sd(get_list_par_from_fit(best_estimation_fits))
  
  list(curve_fits = best_estimation_fits, summary = estimated_summary)
}

calculate_mean_sd <- function(list_of_vectors) {
  if (length(list_of_vectors) == 0) return(list(mean = numeric(0), sd = numeric(0)))
  
  values_matrix <- do.call(rbind, list_of_vectors)
  
  list(
    mean = colMeans(values_matrix),
    sd = apply(values_matrix, 2, sd)
  )
}

get_list_par_from_fit <- function(list_of_fit_objects) {
  lapply(list_of_fit_objects, `[[`, "par")
}
