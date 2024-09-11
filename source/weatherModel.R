# R source code for the demonstration of
# the Indus Village model - Weather model

# visit this model repository at https://github.com/Andros-Spica/indus-village-model/blob/master/01-weather

#=======================================================================
 
# base generic functions

## on solar radiation and temperature

################################################################################
#' @title Calculate the value for a given day in an annual sinusoidal curve
#' @description 
#' Calculate the value for a given day in an annual sinusoidal curve 
#' defined by maximum and minimum values, the length of the year in days, 
#' and the the day with the lowest value in an year.
#' 
#' @param day_of_year Integer. The day of year for which to calculate the value.
#' @param min_value Numeric. The minimum value of the sinusoidal curve within the year.
#' @param max_value Numeric. The maximum value of the sinusoidal curve within the year.
#' @param year_length Integer. The number of days in the given year.
#' @param lowest_value_day Integer. The day with the lowest value in a year.
#'
#' @return Numeric. Point value for the day of year (same units as min_value and max_value).
#' @export
#'
#' @examples
#' get_day_value_in_annual_sinusoid(180, 0, 100, 365, 1)
get_day_value_in_annual_sinusoid <- function(day_of_year,
                                             min_value, 
                                             max_value,  
                                             year_length, 
                                             lowest_value_day) {
  amplitude <- (max_value - min_value) / 2
  phase_shift <- 2 * pi * (day_of_year - lowest_value_day + year_length) / year_length - pi / 2
  
  min_value + amplitude * (1 + sin(phase_shift))
}

################################################################################
#' @title Get day of year with lowest value in annual sinusoidal curve
#' @description 
#' get the day of year with the lowest value in annual sinusoidal curve 
#' according to whether it refers to Earth's southern hemisphere.
#' 
#' @param is_southern_hemisphere Logical. TRUE if the annual curve corresponds to 
#'   values in the southern hemisphere, FALSE for the northern hemisphere.
#' @param year_length Numeric. The number of days in the given year.
#'
#' @return Integer. Day of year with the lowest value in the annual sinusoidal curve.
#' @export
#'
#' @examples
#' get_lowest_value_day(is_southern_hemisphere = FALSE, year_length = 365)
#' get_lowest_value_day(is_southern_hemisphere = TRUE, year_length = 366)
get_lowest_value_day <- function(is_southern_hemisphere, year_length) {
  if (is_southern_hemisphere) {
    lowest_day_fraction <- 0.471  # Approximately June 21st (172 / 365 = 0471...)
  } else {
    lowest_day_fraction <- 0.972  # Approximately December 21st (355 / 365 = 0.972...)
  }
  
  round(lowest_day_fraction * year_length)
}

################################################################################
#' @title Get a day's value in an annual sinusoidal curve with fluctuations
#' @description 
#' Calculate the value for a given day in an annual sinusoidal curve 
#' defined by maximum and minimum values, the length of the year in days, 
#' and whether it refers to Earth's southern hemisphere.
#' 
#' @param day_of_year Integer. The day of year for which to calculate the value.
#' @param min_value Numeric. Minimum value of the sinusoidal curve within the year.
#' @param max_value Numeric. Maximum value of the sinusoidal curve within the year.
#' @param year_length Integer. The number of days in the given year.
#' @param is_southern_hemisphere Logical. TRUE if the curve corresponds to the 
#'   southern hemisphere, FALSE for northern. Default is FALSE.
#' @param fluctuation Numeric. Standard deviation of the normal random noise to be 
#'   added to the sinusoidal curve (same units as min_value and max_value).
#' @param seed Integer. Random number generator seed. Default is NULL, assuming seed is already set.
#'
#' @return Numeric. Point value for the day of year (same units as min_value and max_value).
#' @export
#'
#' @examples
#' get_day_annual_sinusoid_with_fluctuation(180, 0, 100, 365, FALSE, 5, 123)
get_day_annual_sinusoid_with_fluctuation <- function(day_of_year, 
                                                     min_value, 
                                                     max_value, 
                                                     year_length, 
                                                     is_southern_hemisphere = FALSE,
                                                     fluctuation, 
                                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  lowest_value_day <- get_lowest_value_day(is_southern_hemisphere, year_length)
  
  base_value <- get_day_value_in_annual_sinusoid(
    day_of_year,
    min_value, 
    max_value,
    year_length, 
    lowest_value_day
  )
  
  rnorm(1, mean = base_value, sd = fluctuation)
}

################################################################################
#' @title Get an annual sinusoidal curve
#' @description 
#' Calculate an annual sinusoidal curve defined by maximum and minimum values, 
#' the length of the year in days, and the hemisphere (northern or southern).
#' 
#' @param min_value Numeric. Minimum value of the sinusoidal curve within the year.
#' @param max_value Numeric. Maximum value of the sinusoidal curve within the year.
#' @param year_length Integer. The number of days in the given year.
#' @param is_southern_hemisphere Logical. TRUE if the curve corresponds to the 
#'   southern hemisphere, FALSE for northern. Default is FALSE.
#'
#' @return Numeric vector. Values for each day of year describing an annual 
#'   sinusoidal curve (same units as min_value and max_value).
#' @export
#'
#' @examples
#' annual_curve <- get_annual_sinusoid(0, 100, 365, FALSE)
#' plot(annual_curve, type = "l", xlab = "Day of Year", ylab = "Value")
get_annual_sinusoid <- function(min_value, 
                                max_value,  
                                year_length, 
                                is_southern_hemisphere = FALSE) {
  lowest_value_day <- get_lowest_value_day(is_southern_hemisphere, year_length)
  
  sapply(1:year_length, function(day) {
    get_day_value_in_annual_sinusoid(
      day_of_year = day,
      min_value = min_value, 
      max_value = max_value,
      year_length = year_length, 
      lowest_value_day = lowest_value_day
    )
  })
}

################################################################################
#' @title Get an annual sinusoidal curve with fluctuation
#' @description 
#' Calculate an annual sinusoidal curve with random fluctuations, defined by 
#' maximum and minimum values, the length of the year in days, and the hemisphere.
#'
#' @param min_value Numeric. Minimum value of the sinusoidal curve within the year.
#' @param max_value Numeric. Maximum value of the sinusoidal curve within the year.
#' @param year_length Integer. The number of days in the given year.
#' @param is_southern_hemisphere Logical. TRUE if the curve corresponds to the 
#'   southern hemisphere, FALSE for northern. Default is FALSE.
#' @param fluctuation Numeric. Standard deviation of the normal random noise to be 
#'   added to the sinusoidal curve (same units as min_value and max_value).
#' @param seed Integer. Seed for random number generation. Default is NULL, assuming seed is already set.
#'
#' @return Numeric vector. Values for each day of year describing an annual 
#'   sinusoidal curve with fluctuations (same units as min_value and max_value).
#' @export
#'
#' @examples
#' set.seed(123)
#' annual_curve <- get_annual_sinusoid_with_fluctuation(0, 100, 365, FALSE, 5)
#' plot(annual_curve, type = "l", xlab = "Day of Year", ylab = "Value")
get_annual_sinusoid_with_fluctuation <- function(min_value, 
                                                 max_value,  
                                                 year_length, 
                                                 is_southern_hemisphere = FALSE,
                                                 fluctuation,
                                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sapply(1:year_length, function(day) {
    get_day_annual_sinusoid_with_fluctuation(
      day_of_year = day,
      min_value = min_value, 
      max_value = max_value,
      year_length = year_length, 
      is_southern_hemisphere = is_southern_hemisphere,
      fluctuation = fluctuation
    )
  })
}

#=======================================================================

## on precipitation

################################################################################
#' @title Get a day's value in an annual double logistic curve
#' @description 
#' Calculate the value for a given day in an annual double logistic curve 
#' defined by five shape parameters.
#' 
#' @param day_of_year Numeric. The day of the year (1-365 or 1-366).
#' @param plateau_value Numeric. The value (range of 0 to 1) at which the gap between logistic curves is set.
#' @param inflection1 Numeric. The day of year when the first logistic curve has its maximum slope.
#' @param rate1 Numeric. The maximum rate or slope increase of the first logistic curve.
#' @param inflection2 Numeric. The day of year when the second logistic curve has its maximum slope.
#' @param rate2 Numeric. The maximum rate or slope increase of the second logistic curve.
#'
#' @return Numeric. Point value for the day of year (range of 0 to 1).
#' @export
#'
#' @examples
#' get_day_value_in_annual_double_logistic(180, 0.5, 100, 0.1, 250, 0.1)
get_day_value_in_annual_double_logistic <- function(day_of_year, 
                                                    plateau_value, 
                                                    inflection1, 
                                                    rate1, 
                                                    inflection2, 
                                                    rate2) {
  validate_inputs(day_of_year, plateau_value, inflection1, rate1, inflection2, rate2)
  
  logistic1 <- calculate_logistic(day_of_year, inflection1, rate1)
  logistic2 <- calculate_logistic(day_of_year, inflection2, rate2)
  
  plateau_value * logistic1 + (1 - plateau_value) * logistic2
}

#' @title Calculate a single logistic curve value
#'
#' @param x Numeric. The input value.
#' @param inflection Numeric. The inflection point of the logistic curve.
#' @param rate Numeric. The rate of change of the logistic curve.
#'
#' @return Numeric. The value of the logistic curve at point x.
calculate_logistic <- function(x, inflection, rate) {
  1 / (1 + exp((inflection - x) * rate))
}

#' @title Validate input parameters
#'
#' @param day_of_year Numeric. The day of the year (1-365 or 1-366).
#' @param plateau_value Numeric. The value (range of 0 to 1) at which the gap between logistic curves is set.
#' @param inflection1,inflection2 Numeric. The days of year when the logistic curves have their maximum slope.
#' @param rate1,rate2 Numeric. The maximum rates or slope increases of the logistic curves.
#'
#' @return NULL. Throws an error if inputs are invalid.
validate_inputs <- function(day_of_year, plateau_value, inflection1, rate1, inflection2, rate2) {
  if (!is.numeric(day_of_year) || day_of_year < 1 || day_of_year > 366) {
    stop("day_of_year must be a number between 1 and 366")
  }
  if (!is.numeric(plateau_value) || plateau_value < 0 || plateau_value > 1) {
    stop("plateau_value must be a number between 0 and 1")
  }
  if (!is.numeric(inflection1) || !is.numeric(inflection2) || 
      inflection1 < 1 || inflection1 > 366 || inflection2 < 1 || inflection2 > 366) {
    stop("inflection1 and inflection2 must be numbers between 1 and 366")
  }
  if (!is.numeric(rate1) || !is.numeric(rate2)) {
    stop("rate1 and rate2 must be numeric")
  }
}

################################################################################
#' @title Calculate the daily values of an annual double logistic curve
#' @description 
#' Calculate the daily values of an annual double logistic curve defined 
#' by five shape parameters.
#' 
#' @param plateau_value Numeric. The value (range of 0 to 1) at which the gap between logistic curves is set.
#' @param inflection1 Numeric. The day of year when the first logistic curve has its maximum slope.
#' @param rate1 Numeric. The maximum rate or slope increase of the first logistic curve.
#' @param inflection2 Numeric. The day of year when the second logistic curve has its maximum slope.
#' @param rate2 Numeric. The maximum rate or slope increase of the second logistic curve.
#' @param year_length Numeric. The number of days in the given year.
#'
#' @return Numeric vector. Daily values following an annual double logistic curve (range of 0 to 1).
#' @export
#'
#' @examples
#' curve <- get_annual_double_logistic_curve(0.5, 100, 0.1, 250, 0.1, 365)
get_annual_double_logistic_curve <- function(plateau_value,
                                             inflection1,
                                             rate1,
                                             inflection2,
                                             rate2,
                                             year_length) {
  validate_inputs(year_length, plateau_value, inflection1, rate1, inflection2, rate2)
  
  days <- seq_len(year_length)
  vapply(days, get_day_value_in_annual_double_logistic, 
         FUN.VALUE = numeric(1),
         plateau_value = plateau_value, 
         inflection1 = inflection1, 
         rate1 = rate1, 
         inflection2 = inflection2, 
         rate2 = rate2)
}

################################################################################
#' @title Discretise curve stochastically
#'
#' @description Break curve slope into several random steps, each consisting of an increase with maximum slope followed by a plateau.
#'
#' @param curve Numeric vector. The y variable describing the curve to be modified.
#' @param n_samples Integer. The number of random samples or steps to be created in the curve.
#' @param max_sample_size Integer. The maximum length of samples or step plateaus.
#' @param seed Integer. Random number generator seed. Default is NULL, assuming seed is already set.
#'
#' @return Numeric vector. A discretised version of the curve given as input.
#' @export
#'
#' @examples
#' set.seed(123)
#' original_curve <- runif(100)
#' discretised_curve <- discretise_curve(original_curve, n_samples = 10, max_sample_size = 5)
discretise_curve <- function(curve, n_samples, max_sample_size, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  curve_length <- length(curve)
  
  for (i in seq_len(n_samples)) {
    sample_size <- ceiling(max_sample_size * i / n_samples)
    plateau_middle_point <- sample(curve_length, 1)
    
    earliest_neighbor <- max(1, plateau_middle_point - sample_size)
    latest_neighbor <- min(curve_length, plateau_middle_point + sample_size)
    
    neighborhood <- earliest_neighbor:latest_neighbor
    mean_neighborhood <- mean(curve[neighborhood])
    
    curve[neighborhood] <- mean_neighborhood
  }
  
  curve
}

################################################################################
#' @title Rescale curve to 0-1 interval
#' @description Rescale curve to 0-1 interval. If curve is flat, returns original curve.
#' @param curve Numeric vector. The y variable describing the curve to be modified.
#'
#' @return Numeric vector. The rescaled curve in the 0-1 interval.
#' @export
#'
#' @examples
#' original_curve <- c(1, 2, 3, 4, 5)
#' rescaled_curve <- rescale_curve(original_curve)
rescale_curve <- function(curve) {
  if (length(curve) < 2) {
    stop("Curve must have at least two points")
  }
  
  range <- curve[length(curve)] - curve[1]
  
  if (range == 0) {
    warning("Curve has constant value, returning original curve")
    return(curve)
  }
  
  (curve - curve[1]) / range
}

################################################################################
#' @title Get increments from cumulative curve
#' @description Calculate the incremental or differences curve corresponding to a given cumulative curve
#' @param cumulative_curve Numeric vector. The cumulative curve from which to calculate the incremental curve.
#'
#' @return Numeric vector. The incremental curve corresponding to the given cumulative curve.
#' @export
#'
#' @examples
#' cumulative_curve <- c(1, 3, 6, 10, 15)
#' incremental_curve <- get_increments_from_cumulative_curve(cumulative_curve)
get_increments_from_cumulative_curve <- function(cumulative_curve) {
  if (length(cumulative_curve) < 2) {
    stop("Cumulative curve must have at least two points")
  }
  
  increment_curve <- c(cumulative_curve[1], diff(cumulative_curve))
  
  # Correct negative values (in rare occasions, very small negative numbers may appear)
  pmax(increment_curve, 0)
}

################################################################################
#' @title Get daily precipitation of one year
#'
#' @description Calculate the daily values of precipitation using the double logistic cumulative curve approach.
#'
#' @param plateau_value Numeric. The value (range of 0 to 1) in which the gap between logistic curves is set.
#' @param inflection1,inflection2 Numeric. The days of year in which the first and second logistic curves have their maximum slope.
#' @param rate1,rate2 Numeric. The maximum rate or slope increase of the first and second logistic curves.
#' @param year_length Integer. The number of days in the given year.
#' @param n_samples Integer. The number of random samples or steps to be created in the cumulative curve.
#' @param max_sample_size Integer. The maximum length of samples or step plateaus in the cumulative curve.
#' @param annual_sum Numeric. The annual sum of precipitation (mm).
#' @param seed Integer. Random number generator seed. Default is NULL, assuming seed is already set.
#'
#' @return Numeric vector. Daily precipitation values (mm).
#' @export
#'
#' @examples
#' precipitation <- get_precipitation_of_year(
#'   plateau_value = 0.5, inflection1 = 100, rate1 = 0.05,
#'   inflection2 = 250, rate2 = 0.05, year_length = 365,
#'   n_samples = 10, max_sample_size = 5, annual_sum = 1000, seed = 123
#' )
get_precipitation_of_year <- function(plateau_value, 
                                      inflection1, 
                                      rate1, 
                                      inflection2, 
                                      rate2, 
                                      year_length, 
                                      n_samples, 
                                      max_sample_size, 
                                      annual_sum,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  precipitation <- get_annual_double_logistic_curve(
    plateau_value = plateau_value,
    inflection1 = inflection1,
    rate1 = rate1,
    inflection2 = inflection2,
    rate2 = rate2,
    year_length = year_length
  )
  
  precipitation <- discretise_curve(
    curve = precipitation,
    n_samples = n_samples, 
    max_sample_size = max_sample_size
  )
  
  # Handle special case where the curve is a horizontal line (first = last)
  if (precipitation[1] == precipitation[length(precipitation)]) { 
    precipitation <- seq_len(year_length) / year_length
  }
  
  # Check if discretisation was sufficient
  if (utils::tail(precipitation, 1) < 1) { 
    warning("Failed to generate precipitation that fulfills 'annual_sum' without re-scaling")
  }
  
  precipitation <- rescale_curve(precipitation)
  precipitation <- get_increments_from_cumulative_curve(precipitation) * annual_sum
  
  precipitation
}

################################################################################
#' @title Get daily proportion to annual sum of precipitation of one year
#'
#' @param daily_precipitation_year Numeric vector. Daily values of precipitation for a length of one year.
#' @return Numeric vector. Daily proportion of annual precipitation (mm/mm).
#' @export
#'
#' @examples
#' daily_precip <- c(1, 2, 3, 4, 5)
#' get_proportion_precipitation_of_year(daily_precip)
get_proportion_precipitation_of_year <- function(daily_precipitation_year) {
  daily_precipitation_year / sum(daily_precipitation_year)
}

################################################################################
#' @title Get daily proportion to annual sum of precipitation of multiple years
#'
#' @param daily_precipitation Numeric vector. Daily values of precipitation for multiple full years.
#' @param years Integer vector. Year for each daily value in daily_precipitation.
#' @return Numeric vector. Daily proportion of annual precipitation (mm/mm).
#' @export
#'
#' @examples
#' daily_precip <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
#' years <- c(rep(2020, 5), rep(2021, 5))
#' get_proportion_precipitation(daily_precip, years)
get_proportion_precipitation <- function(daily_precipitation, years) {
  vapply(split(daily_precipitation, years), 
         get_proportion_precipitation_of_year, 
         FUN.VALUE = numeric(length(daily_precipitation) / length(unique(years))))
}

################################################################################
#' @title Get daily cumulative sum of precipitation of one year
#'
#' @param daily_precipitation_year Numeric vector. Daily values of precipitation for a length of one year.
#' @return Numeric vector. Daily cumulative sum of annual precipitation (mm/mm).
#' @export
#'
#' @examples
#' daily_precip <- c(1, 2, 3, 4, 5)
#' get_cumulative_precipitation_of_year(daily_precip)
get_cumulative_precipitation_of_year <- function(daily_precipitation_year) {
  cumsum(daily_precipitation_year) / sum(daily_precipitation_year)
}

################################################################################
#' @title Get daily cumulative sum of precipitation of multiple years
#'
#' @param daily_precipitation Numeric vector. Daily values of precipitation for multiple full years.
#' @param years Integer vector. Year for each daily value in daily_precipitation.
#' @return Numeric vector. Daily cumulative sum of annual precipitation (mm/mm).
#' @export
#'
#' @examples
#' daily_precip <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
#' years <- c(rep(2020, 5), rep(2021, 5))
#' get_cumulative_precipitation(daily_precip, years)
get_cumulative_precipitation <- function(daily_precipitation, years) {
  unlist(lapply(split(daily_precipitation, years), 
                get_cumulative_precipitation_of_year))
}

#=======================================================================

## reference evapotranspiration
################################################################################
#' @title Estimate daily reference evapotranspiration (ETr)
#'
#' @description 
#' Estimates daily reference evapotranspiration using either Penman-Monteith equation (default) or Priestley-Taylor equation.
#'
#' @param R_s Numeric. Solar radiation or insolation (MJ m^2 day^-1)
#' @param temperature Numeric. Daily average temperature at 2m (°C)
#' @param temperature_max Numeric. Daily maximum temperature at 2m (°C)
#' @param temperature_min Numeric. Daily minimum temperature at 2m (°C)
#' @param temperature_dew Numeric. Daily dew point temperature at 2m (°C). Default is temperature_min.
#' @param wind_speed Numeric. Mean daily wind speed at 2m (m s^-1). Default is 2.
#' @param albedo Numeric. Canopy reflection or albedo. Default is 0.23.
#' @param z Numeric. Elevation above sea level (m). Default is 200.
#' @param lambda Numeric. Latent heat of vaporisation (MJ kg^-1). Default is 2.45.
#' @param c_p Numeric. Specific heat at constant pressure (MJ kg^-1 °C^-1). Default is 1.013e-3.
#' @param epsilon Numeric. Ratio of molecular weight of water vapour to dry air. Default is 0.622.
#' @param method Character. "PM" for Penman-Monteith equation, "PT" for Priestley-Taylor equation. Default is "PM".
#' @param C_n Numeric. Constant used in Penman-Monteith equation. Default is 900 (grass reference).
#' @param C_d Numeric. Constant used in Penman-Monteith equation. Default is 0.34 (grass reference).
#' @param alpha Numeric. Priestley-Taylor coefficient. Default is 1.26.
#'
#' @details # useful references:
#' Suleiman A A and Hoogenboom G 2007 
#' Comparison of Priestley-Taylor and FAO-56 Penman-Monteith for Daily Reference Evapotranspiration Estimation in Georgia 
#' J. Irrig. Drain. Eng. 133 175–82 Online: http://ascelibrary.org/doi/10.1061/%28ASCE%290733-9437%282007%29133%3A2%28175%29
#' also: Jia et al. 2013 - doi:10.4172/2168-9768.1000112
#' constants found in: http://www.fao.org/3/X0490E/x0490e07.htm
#' see also r package: Evapotranspiration (consult source code).
#' 
#' @return Numeric. Daily ETr (mm day^-1)
#' @export
#'
#' @examples
#' estimate_etr(R_s = 20, temperature = 25, temperature_max = 30, temperature_min = 20)
estimate_etr <- function(R_s, 
                         temperature, temperature_max, temperature_min,
                         temperature_dew = temperature_min,
                         wind_speed = 2, 
                         albedo = 0.23,
                         z = 200,
                         lambda = 2.45,
                         c_p = 1.013e-3,
                         epsilon = 0.622,
                         method = 'PM',
                         C_n = 900,
                         C_d = 0.34,
                         alpha = 1.26) {
  
  # Net solar radiation
  net_solar_radiation <- (1 - albedo) * R_s
  
  # Saturated vapor pressure function
  e_o <- function(temp) 0.6108 * exp(17.27 * temp / (temp + 237.3))
  
  # Saturated and actual vapor pressure
  e_s <- (e_o(temperature_max) + e_o(temperature_min)) / 2
  e_a <- e_o(temperature_dew)
  
  # Slope of the vapor pressure-temperature curve (kPa °C^-1)
  delta <- 4098 * e_o(temperature) / (temperature + 237.3)^2
  
  # Atmospheric pressure (kPa)
  P <- 101.3 * ((293 - 0.0065 * z) / 293)^5.26
  
  # Psychrometric constant (kPa °C^-1)
  gamma <- c_p * P / (epsilon * lambda) 
  
  # Calculate ETr based on method
  ETr <- if (method == 'PM') {
    (0.408 * delta * net_solar_radiation + gamma * (C_n / (temperature + 273)) * wind_speed * (e_s - e_a)) / 
      (delta + gamma * (1 + C_d * wind_speed))
  } else if (method == 'PT') {
    (alpha / lambda) * (delta / (delta + gamma)) * net_solar_radiation
  } else {
    stop("Invalid method. Use 'PM' for Penman-Monteith or 'PT' for Priestley-Taylor.")
  }
  
  ETr
}

#=======================================================================

# auxiliar functions

################################################################################
#' @title Get last item in vector
#'
#' @param x A vector of any type
#' @return The last item in the vector
#' @export
#'
#' @examples
#' get_last_item(c(1, 2, 3, 4, 5))
get_last_item <- function(x) {
  if (length(x) == 0) {
    stop("Input vector is empty")
  }
  x[length(x)]
}

################################################################################
#' @title Create vertical lines marking the end of years in a time-series plot
#' @description Calculate the incremental or differences curve corresponding to a given cumulative curve
#' @param length_of_data Integer. Length of the x-axis in the data
#' @param offset Numeric. Scaling offset of each mark (helpful to adjust marks if barplot is used). Default is 1.
#' @param year_length Integer. The number of days per year. Default is 365.
#' @param lty Integer. Line type for vertical lines. Default is 3 (dotted).
#'
#' @return NULL. This function is called for its side effect of adding vertical lines to an existing plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plot(1:1000, rnorm(1000), type = "l")
#' mark_end_years(1000)
#' }
mark_end_years <- function(length_of_data, 
                           offset = 1,
                           year_length = 365,
                           lty = 3) {
  if (!is.numeric(length_of_data) || length_of_data <= 0) {
    stop("length_of_data must be a positive number")
  }
  
  year_ends <- seq(year_length * offset, length_of_data, by = year_length * offset)
  abline(v = year_ends, lty = lty)
}

#=======================================================================

# Main functions

################################################################################
#' @title Initialize the Weather Model
#'
#' @description 
#' Creates a complex R object containing parameter values and empty holders for variables 
#' to initialize the Weather model.
#'
#' @param year_length Integer. Number of days per year. Default is 365.
#' @param seed Integer. Seed for the random number generator. Default is 0.
#' @param albedo Numeric. Albedo used in the calculation of reference evapotranspiration (range 0 to 1). Default is 0.4.
#' @param is_southern_hemisphere Logical. Whether the annual sinusoidal curves correspond to the southern hemisphere. Default is FALSE.
#' @param temp_annual_max Numeric. Annual maximum daily average temperature at 2m (°C). Default is 40.
#' @param temp_annual_min Numeric. Annual minimum daily average temperature at 2m (°C). Default is 15.
#' @param temp_daily_fluctuation Numeric. Standard deviation of daily temperature fluctuation (°C). Default is 5.
#' @param temp_daily_lower_dev Numeric. Daily minimum temperature deviation (°C). Default is 5.
#' @param temp_daily_upper_dev Numeric. Daily maximum temperature deviation (°C). Default is 5.
#' @param solar_annual_max Numeric. Annual maximum daily solar radiation (MJ/m^2). Default is 7.
#' @param solar_annual_min Numeric. Annual minimum daily solar radiation (MJ/m^2). Default is 3.
#' @param solar_daily_fluctuation Numeric. Standard deviation of daily solar radiation fluctuation (MJ/m^2). Default is 1.
#' @param precip_yearly_mean Numeric. Annual sum of precipitation (mm). Default is 400.
#' @param precip_yearly_sd Numeric. Standard deviation of annual precipitation (mm). Default is 130.
#' @param precip_params Named list. Contains parameters for precipitation model. See details in function body.
#'
#' @return A list containing three main components: PARS (parameter values), 
#' annualPrecipitationPars (parameters for annual precipitation), and daily (empty vectors for time-series variables).
#' 
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model(seed = 123, temp_annual_max = 35, precip_yearly_mean = 500)
initialize_weather_model <- function(
  year_length = 365,
  seed = 0,
  albedo = 0.4,
  is_southern_hemisphere = FALSE,
  temp_annual_max = 40,
  temp_annual_min = 15,
  temp_daily_fluctuation = 5,
  temp_daily_lower_dev = 5,
  temp_daily_upper_dev = 5,
  solar_annual_max = 7,
  solar_annual_min = 3,
  solar_daily_fluctuation = 1,
  precip_yearly_mean = 400,
  precip_yearly_sd = 130,
  precip_params = list(
    n_samples_mean = 200,
    n_samples_sd = 5,
    max_sample_size_mean = 10,
    max_sample_size_sd = 3,
    plateau_value_mean = 0.1,
    plateau_value_sd = 0.05,
    inflection1_mean = 40,
    inflection1_sd = 20,
    rate1_mean = 0.15,
    rate1_sd = 0.02,
    inflection2_mean = 200,
    inflection2_sd = 20,
    rate2_mean = 0.05,
    rate2_sd = 0.01
  )
) {
  set.seed(seed)
  
  list(
    PARS = list(
      seed = seed,
      year_length = year_length,
      albedo = albedo,
      is_southern_hemisphere = is_southern_hemisphere,
      temperature = list(
        annual_max = temp_annual_max,
        annual_min = temp_annual_min,
        daily_fluctuation = temp_daily_fluctuation,
        daily_lower_dev = temp_daily_lower_dev,
        daily_upper_dev = temp_daily_upper_dev
      ),
      solar = list(
        annual_max = solar_annual_max,
        annual_min = solar_annual_min,
        daily_fluctuation = solar_daily_fluctuation
      ),
      precipitation = c(
        list(yearly_mean = precip_yearly_mean, yearly_sd = precip_yearly_sd),
        precip_params
      )
    ),
    annual_precipitation_pars = list(
      annual_sum = numeric(0),
      plateau_value = numeric(0),
      inflection1 = numeric(0),
      rate1 = numeric(0),
      inflection2 = numeric(0),
      rate2 = numeric(0),
      n_samples = numeric(0),
      max_sample_size = numeric(0)
    ),
    daily = list(
      current_year = integer(0),
      current_day_of_year = integer(0),
      temperature = numeric(0),
      temperature_max = numeric(0),
      temperature_min = numeric(0),
      solar_radiation = numeric(0),
      ETr = numeric(0),
      precipitation = numeric(0)
    )
  )
}

################################################################################
#' @title Update Precipitation Parameters
#'
#' @description 
#' Updates the precipitation parameters for the given weather model instance
#'
#' @param weather_model List. Initialized instance of weather model
#' @param year_length Integer. Number of days in a year
#'
#' @return Updated weather model with new precipitation parameters
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model()
#' weather_model <- update_precipitation_parameters(weather_model, 365)
update_precipitation_parameters <- function(weather_model, year_length) {
  precip_pars <- weather_model$PARS$precipitation
  
  new_params <- list(
    annual_sum = max(0, rnorm(1, precip_pars$yearly_mean, precip_pars$yearly_sd)),
    plateau_value = pmin(1, pmax(0, rnorm(1, precip_pars$plateau_value_mean, precip_pars$plateau_value_sd))),
    inflection1 = pmin(year_length, pmax(1, rnorm(1, precip_pars$inflection1_mean, precip_pars$inflection1_sd))),
    rate1 = max(0, rnorm(1, precip_pars$rate1_mean, precip_pars$rate1_sd)),
    inflection2 = pmin(year_length, pmax(1, rnorm(1, precip_pars$inflection2_mean, precip_pars$inflection2_sd))),
    rate2 = max(0, rnorm(1, precip_pars$rate2_mean, precip_pars$rate2_sd)),
    n_samples = max(0, rnorm(1, precip_pars$n_samples_mean, precip_pars$n_samples_sd)),
    max_sample_size = max(0, rnorm(1, precip_pars$max_sample_size_mean, precip_pars$max_sample_size_sd))
  )
  
  for (param_name in names(new_params)) {
    weather_model$annual_precipitation_pars[[param_name]] <- c(
      weather_model$annual_precipitation_pars[[param_name]],
      new_params[[param_name]]
    )
  }
  
  weather_model
}

################################################################################
#' @title Generate Annual Precipitation
#'
#' @description 
#' Generates annual precipitation values based on the updated parameters in the given weather model instance
#'
#' @param weather_model List. Initialized instance of weather model
#' @param year_length Integer. Number of days in a year
#'
#' @return Updated weather model with new precipitation values
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model()
#' weather_model <- update_precipitation_parameters(weather_model, 365)
#' weather_model <- generate_annual_precipitation(weather_model, 365)
generate_annual_precipitation <- function(weather_model, year_length) {
  last_params <- lapply(weather_model$annual_precipitation_pars, tail, n = 1)
  
  new_precipitation <- get_precipitation_of_year(
    plateau_value = last_params$plateau_value,
    inflection1 = last_params$inflection1,
    rate1 = last_params$rate1,
    inflection2 = last_params$inflection2,
    rate2 = last_params$rate2,
    year_length = year_length,
    n_samples = last_params$n_samples,
    max_sample_size = last_params$max_sample_size,
    annual_sum = last_params$annual_sum
  )
  
  weather_model$daily$precipitation <- c(weather_model$daily$precipitation, new_precipitation)
  
  weather_model
}

################################################################################
#' @title Generate Temperature Variables for a given day
#'
#' @description 
#' Generates mean, min, and maximum temperature for a given day in the year in the given weather model instance
#'
#' @param weather_model List. Initialized instance of weather model
#' @param day Integer. Day of year
#' @param year_length Integer. Number of days in the year
#'
#' @return Updated weather model with new temperature variables
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model()
#' weather_model <- generate_temperature_variables(weather_model, 1, 365)
generate_temperature_variables <- function(weather_model, day, year_length) {
  temp_params <- weather_model$PARS$temperature
  
  mean_temp <- get_day_annual_sinusoid_with_fluctuation(
    min_value = temp_params$annual_min,
    max_value = temp_params$annual_max,
    fluctuation = temp_params$daily_fluctuation,
    day_of_year = day,
    year_length = year_length,
    is_southern_hemisphere = weather_model$PARS$is_southern_hemisphere
  )
  
  weather_model$daily$temperature <- c(weather_model$daily$temperature, mean_temp)
  weather_model$daily$temperature_min <- c(weather_model$daily$temperature_min, mean_temp - temp_params$daily_lower_dev)
  weather_model$daily$temperature_max <- c(weather_model$daily$temperature_max, mean_temp + temp_params$daily_upper_dev)
  
  weather_model
}

################################################################################
#' @title Generate Solar Radiation for a given day
#'
#' @description 
#' Generates total solar radiation for a given day in the year in the given weather model instance
#'
#' @param weather_model List. Initialized instance of weather model
#' @param day Integer. Day of year
#' @param year_length Integer. Number of days in the year
#'
#' @return Updated weather model with new solar radiation value
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model()
#' weather_model <- generate_solar_radiation(weather_model, 1, 365)
generate_solar_radiation <- function(weather_model, day, year_length) {
  solar_params <- weather_model$PARS$solar
  
  solar_radiation <- max(0, get_day_annual_sinusoid_with_fluctuation(
    min_value = solar_params$annual_min,
    max_value = solar_params$annual_max,
    fluctuation = solar_params$daily_fluctuation,
    day_of_year = day,
    year_length = year_length,
    is_southern_hemisphere = weather_model$PARS$is_southern_hemisphere
  ))
  
  weather_model$daily$solar_radiation <- c(weather_model$daily$solar_radiation, solar_radiation)
  
  weather_model
}

################################################################################
#' @title Generate Reference Evapotranspiration for a given day
#'
#' @description 
#' Generates reference evapotranspiration for a given day in the year in the given weather model instance
#'
#' @param weather_model List. Initialized instance of weather model
#'
#' @return Updated weather model with new evapotranspiration value
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model()
#' weather_model <- generate_etr(weather_model)
generate_etr <- function(weather_model) {
  daily <- weather_model$daily
  etr <- estimate_etr(
    R_s = tail(daily$solar_radiation, 1),
    temperature = tail(daily$temperature, 1),
    temperature_max = tail(daily$temperature_max, 1),
    temperature_min = tail(daily$temperature_min, 1)
  )
  
  weather_model$daily$ETr <- c(daily$ETr, etr)
  
  weather_model
}

################################################################################
#' Run the Weather Model
#'
#' @description 
#' Runs the weather model for the given number of years
#'
#' @param weather_model List. Initialized instance of weather model
#' @param num_years Integer. Number of years to be simulated
#'
#' @return Updated weather model containing the corresponding time-series recorded inside the "daily" partition
#' @export
#'
#' @examples
#' weather_model <- initialize_weather_model()
#' weather_model <- run_weather_model(weather_model, 10)
run_weather_model <- function(weather_model, num_years) {
  for (year in seq_len(num_years)) {
    year_length <- weather_model$PARS$year_length[min(year, length(weather_model$PARS$year_length))]
    
    weather_model <- update_precipitation_parameters(weather_model, year_length)
    weather_model <- generate_annual_precipitation(weather_model, year_length)
    
    for (day in seq_len(year_length)) {
      weather_model$daily$current_year <- c(weather_model$daily$current_year, year)
      weather_model$daily$current_day_of_year <- c(weather_model$daily$current_day_of_year, day)
      
      weather_model <- generate_temperature_variables(weather_model, day, year_length)
      weather_model <- generate_solar_radiation(weather_model, day, year_length)
      weather_model <- generate_etr(weather_model)
    }
  }
  
  weather_model
}
