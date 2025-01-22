print_parameter_comparison_table <- function(original_params, fit, params_range_upper, params_range_lower) {
  # Create parameter comparison data frame
  parameter_comparison <- cbind(
    original = original_params,
    estimated = fit$par,
    delta = round(abs(original_params - fit$par), digits = 6),
    range = params_range_upper - params_range_lower,
    delta_percent = round(
      100 * abs(original_params - fit$par) / (params_range_upper - params_range_lower), 
      digits = 4
    )
  )
  
  knitr::kable(parameter_comparison, 
               format = "html",
               col.names = c("original", "estimated", "delta", 
                             "range", "delta (%)"),
               align = c("c", "c", "c", "c", "c"))
}
