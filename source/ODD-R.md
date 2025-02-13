# ODD document for the Indus Village's Weather model (R implementation)

## Purpose and Patterns

The Weather model, implemented in R, is designed to generate daily weather variables for use in simulation models. It simulates weather patterns, including daily temperature, precipitation, and solar radiation. The primary goal is to provide realistic weather data that can drive other simulations, such as those modeling ecological or agricultural systems. The R code includes functions for generating sinusoidal curves for temperature and solar radiation, as well as more complex functions for simulating precipitation patterns.

## Entities, State Variables, and Scales

**Entities**: The model primarily deals with environmental variables at a specific location. There are no explicit agents or individual entities in this model.

**State Variables**:

*   `temperature`: Daily average temperature (°C).

*   `temperature_max`: Daily maximum temperature (°C).

*   `temperature_min`: Daily minimum temperature (°C).

*   `solar_radiation`: Daily solar radiation (MJ m^-2).

*   `ETr`: Daily reference evapotranspiration (mm day^-1).

*   `precipitation`: Daily precipitation (mm).

*   `current_year`: The current year of the simulation.

*   `current_day_of_year`: The current day of the year.

**Scales**:

*   **Temporal**: The model operates on a daily time step, generating weather variables for each day. It can simulate weather patterns over multiple years. The `year_length` parameter determines the number of days in a year.

*   **Spatial**: The model assumes local, homogeneous conditions. The generated weather variables apply uniformly across the simulated area.

## Process Overview and Scheduling

The R implementation of the Weather model follows a clear sequence of steps to generate daily weather variables:

1.  **Initialization (`initialise_weather_model`)**: The `initialise_weather_model` function sets up the model by defining parameter values and creating empty vectors to store daily weather variables. Key parameters include annual maximum and minimum temperatures (`temp_annual_max`, `temp_annual_min`), solar radiation parameters (`solar_annual_max`, `solar_annual_min`), and precipitation parameters (`precip_annual_sum_mean`, `precip_annual_sum_sd`). The function returns a list containing parameter values (`PARAMS`), annual precipitation parameters (`annual_precipitation_params`), and empty vectors for daily time-series variables (`daily`).

2.  **Update Precipitation Parameters (`update_precipitation_parameters`)**: This function updates precipitation parameters for each year. It introduces variability by sampling from normal distributions, ensuring that parameters like `annual_sum`, `plateau_value`, `inflection1`, `rate1`, `inflection2`, `rate2`, `n_samples`, and `max_sample_size` are updated.

3.  **Generate Annual Precipitation (`generate_annual_precipitation`)**: This function generates daily precipitation values for the year based on the updated precipitation parameters. It uses the `gen_precipitation_of_year` function, which employs a double logistic curve to simulate the annual precipitation pattern. The `rescale_curve` function ensures the cumulative precipitation curve is scaled to the 0-1 interval.

4.  **Generate Temperature Variables (`generate_temperature_variables`)**: This function calculates the mean, minimum, and maximum temperatures for a given day. It uses the `gen_day_annual_sinusoid_with_fluctuation` function to generate the mean temperature, incorporating annual sinusoidal patterns and daily fluctuations. Minimum and maximum temperatures are then calculated by subtracting and adding deviations from the mean temperature.

5.  **Generate Solar Radiation (`generate_solar_radiation`)**: This function generates the total solar radiation for a given day. Similar to temperature, it uses the `gen_day_annual_sinusoid_with_fluctuation` function to simulate the annual solar radiation pattern with daily fluctuations.

6.  **Generate Reference Evapotranspiration (`generate_etr`)**: This function estimates the daily reference evapotranspiration (`ETr`) using the `estimate_etr` function. The `estimate_etr` function uses either the Penman-Monteith or Priestley-Taylor equation, incorporating variables such as solar radiation, temperature, wind speed, and albedo.

7.  **Run Weather Model (`run_weather_model`)**: The `run_weather_model` function orchestrates the entire simulation. For each year, it updates precipitation parameters, generates annual precipitation, and then calculates daily temperature, solar radiation, and evapotranspiration. It stores the results in the `daily` list within the weather model object.

## Design Concepts

**Stochasticity**: The R model incorporates stochasticity by using random number generators to introduce variability in precipitation parameters and daily fluctuations in temperature and solar radiation. The `rnorm` function is used to sample from normal distributions.

**Seasonality**: The model captures the seasonal patterns of temperature and solar radiation using sinusoidal functions. The `gen_day_value_in_annual_sinusoid` function calculates the base value for a given day based on the annual cycle.

**Empirical parameters**: The evapotranspiration submodel depends on empirical parameters.

**Emergence**: Daily weather patterns emerge from the combination of sinusoidal functions, stochastic processes, and empirical relationships.

**Sensing**: The model does not include a sensing component.

**Adaptation**: The model does not include an adaptation component.

**Objectives**: The model does not include agent objectives.

**Learning**: The model does not include a learning component.

**Prediction**: The model is primarily a data-generating tool rather than a predictive model.

**Collectives**: The model does not simulate collectives.

**Observation**: The primary output of the model is the set of daily weather variables, which can be used as inputs for other simulation models.

**Initialisation**: The `initialise_weather_model` function initializes the state variables and parameters at the beginning of the simulation.

**Input Data**: The model relies on parameters set at the beginning of the simulation rather than external input data.

## Initialisation

The `initialise_weather_model` function initializes the model. It sets default values for parameters such as `year_length`, `albedo`, `temp_annual_max`, `temp_annual_min`, `solar_annual_max`, and `solar_annual_min`. It also creates empty vectors for storing the daily weather variables.

## Input Data

The model does not rely on external input data. Instead, it uses parametric equations and stochastic processes to generate daily weather variables. The parameters for these equations are set either through default values or by sampling from distributions.

Here is a full account of the constants and modified parameters:

*   `year_length`: Specifies the **number of days per year**; the default is 365.  
*   `seed`: Sets the **seed for the random number generator**; the default is 0.  
*   `albedo`: Defines the **albedo** used to calculate reference evapotranspiration, with a range from 0 to 1; the default value is 0.  4.  
*   `is_southern_hemisphere`: A logical parameter that indicates **whether the annual sinusoidal curves correspond to the southern hemisphere**; the default is FALSE.  
*   `temp_annual_max`: Indicates the **annual maximum daily average temperature at 2m (°C)**; the default is 40.  
*   `temp_annual_min`: Indicates the **annual minimum daily average temperature at 2m (°C)**; the default is 15.  
*   `temp_daily_fluctuation`: Sets the **standard deviation of daily temperature fluctuation (°C)**; the default is 5.  
*   `temp_daily_lower_dev`: Establishes the **daily minimum temperature deviation (°C)**; the default is 5.  
*   `temp_daily_upper_dev`: Establishes the **daily maximum temperature deviation (°C)**; the default is 5.  
*   `solar_annual_max`: Sets the **annual maximum daily solar radiation (MJ/m^2)**; the default is 7.  
*   `solar_annual_min`: Sets the **annual minimum daily solar radiation (MJ/m^2)**; the default is 3.  
*   `solar_daily_fluctuation`: Sets the **standard deviation of daily solar radiation fluctuation (MJ/m^2)**; the default is 1.  
*   `precip_annual_sum_mean`: Indicates the **mean of annual sum of precipitation (mm)**; the default is 400.  
*   `precip_annual_sum_sd`: Indicates the **standard deviation of annual sum of precipitation (mm)**; the default is 130.  
*   `precip_plateau_value_mean`: Sets the **mean plateau value of the annual cumulative precipitation curve**.   This is the value (range of 0 to 1) in which the gap between logistic curves is set; the default is 0.  1.  
*   `precip_plateau_value_sd`: Sets the **standard deviation of plateau value of the annual cumulative precipitation curve**; the default is 0.05.  
*   `precip_inflection1_mean`: Sets the **mean of the day of year in which the first logistic curve has its maximum slope**; the default is 40.  
*   `precip_inflection1_sd`: Sets the **standard deviation of the days of year in which the first logistic curve has its maximum slope**; the default is 20.  
*   `precip_inflection2_mean`: Sets the **mean of the days of year in which the second logistic curve has its maximum slope**; the default is 200.  
*   `precip_inflection2_sd`: Sets the **standard deviation of the days of year in which the second logistic curve has its maximum slope**; the default is 20.  
*   `precip_rate1_mean`: Sets the **mean of the maximum rate or slope increase of the first logistic curve**; the default is 0.15.  
*   `precip_rate1_sd`: Sets the **standard deviation of the maximum rate or slope increase of the first logistic curve**; the default is 0.02.  
*   `precip_rate2_mean`: Sets the **mean of the maximum rate or slope increase of the second logistic curve**; the default is 0.05.  
*   `precip_rate2_sd`: Sets the **standard deviation of the maximum rate or slope increase of the second logistic curve**; the default is 0.01.  
*   `precip_n_samples_mean`: Sets the **mean of the number of random samples or steps to be created in the cumulative curve**; the default is 200.  
*   `precip_n_samples_sd`: Sets the **standard deviation of the number of random samples or steps to be created in the cumulative curve**; the default is 5.  
*   `precip_max_sample_size_mean`: Sets the **mean of the maximum length of samples or step plateaus in the cumulative curve**; the default is 10.  
*   `precip_max_sample_size_sd`: Sets the **standard deviation of the maximum length of samples or step plateaus in the cumulative curve**; the default is 3.  

## Submodels

### Sinusoidal variables (solar radiation, temperature, atmospheric CO2 concentration)

**`gen_day_value_in_annual_sinusoid`**: This function generates a value for a given day in an annual sinusoidal curve, defined by minimum and maximum values.

**`gen_day_annual_sinusoid_with_fluctuation`**: This function generates a value for a given day in an annual sinusoidal curve, with added random fluctuations.

**`generate_temperature_variables`**: Generates mean, min, and maximum temperature for a given day in the year.

**`generate_solar_radiation`**: Generates total solar radiation for a given day in the year.

(atmospheric CO2 concentration is not implemented in R, because its behaviour is similar to solar radiation; `generate_solar_radiation` can be reused with different arguments for this purpose)

### Precipitation (using double logistic equation)

**`rescale_curve`**: This function rescales a curve to the 0-1 interval. It is used to normalize the cumulative precipitation curve.

**`gen_annual_double_logistic_curve`**: This function generates an annual double logistic curve. It is used to simulate the cumulative precipitation pattern over the year.

**`gen_cum_precipitation_of_year`**: This function generates the cumulative precipitation for the year.

**`gen_precipitation_of_year`**: This function generates the daily precipitation values for the year. It calculates the incremental changes from the cumulative precipitation curve.

**`update_precipitation_parameters`**: Updates the precipitation parameters for the given weather model instance.

### Reference Evapotranspiration (FAO standard)

**`estimate_etr`**: This function estimates the daily reference evapotranspiration using either the Penman-Monteith or Priestley-Taylor equation.

**`generate_etr`**: Generates reference evapotranspiration for a given day.
