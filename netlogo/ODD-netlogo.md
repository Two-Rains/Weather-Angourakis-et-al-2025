# ODD document for the Indus Village's Weather model (NetLogo implementation)

## Purpose and Patterns

The Weather model is designed to procedurally generate daily weather variables for use in larger simulation models. The model aims to produce realistic weather patterns, simulating daily temperature, precipitation, CO2 concentration, and solar radiation.

## Entities, State Variables, and Scales

**Entities**: The model primarily deals with environmental variables at a specific location. There are no explicit agents or individual entities in this model.

**State variables** reflect the current conditions of the simulation at any given time step. These variables are updated as the model progresses and are used to track changes in the simulated environment:

-   **Time tracking**:
    -   `currentYear`: The current year of the simulation.  
    -   `currentDayOfYear`: The current day of the year.  
-   **Main variables**:
    -   `T`: The average temperature of the current day in Celsius.  
    -   `T_max`: The maximum temperature of the current day in Celsius.  
    -   `T_min`: The minimum temperature of the current day in Celsius.  
    -   `RAIN`: The amount of precipitation on the current day in millimeters.  
    -   `CO2`: The average CO2 concentration for the current day in parts per million.  
    -   `solarRadiation`: The amount of solar radiation on the current day in MJ m\^-2.  
    -   `netSolarRadiation`: The net solar radiation, accounting for canopy reflection or albedo.  
    -   `precipitation_yearSeries`: A time series of precipitation values for the entire year.  
    -   `precipitation_cumYearSeries`: A time series of cumulative precipitation values for the year.  
    -   `ETr`: The reference evapotranspiration for the current day.  

**Scales:**

-   **Temporal**: The model operates on a daily time step, generating weather variables for each day of the year. It simulates weather patterns over a single year, with the option to extend the simulation over multiple years.

-   **Spatial**: The spatial scale is not explicitly defined, but the model assumes local homogeneous conditions, meaning that the generated weather variables apply uniformly across the simulated area.

**No breeds**: There are no agent or patch variables.

## Process Overview and Scheduling

**Setup:** The `setup` procedure initializes the model. It calls `set-constants` to define constant values such as `yearLengthInDays`. It calls `set-parameters` to set the model parameters based on the experiment type (user-defined or random). It sets `currentDayOfYear` to 1 and calls `update-weather` to generate the initial weather conditions. Finally, it resets the ticks counter.

**Go:** The `go` procedure simulates the weather patterns over time. It calls `update-weather` to generate the daily weather variables. It calls `advance-time` to increment the `currentDayOfYear` and `currentYear` variables. The simulation continues until the ticks counter reaches the `end-simulation-in-tick` value.

**Update-weather:** The `update-weather` procedure updates the daily weather variables. It calls `update-temperature` and `update-precipitation` to generate the daily temperature and precipitation values. It calls `gen-CO2` and `gen-solar-radiation` to generate the daily CO2 concentration and solar radiation values. It calculates the net solar radiation by discounting the solar radiation based on the albedo. It calculates the reference evapotranspiration (`ETr`) using the `get-ETr` procedure.

**Advance-time:** The `advance-time` procedure increments the `currentDayOfYear` variable by 1. If `currentDayOfYear` exceeds `yearLengthInDays`, it increments the `currentYear` variable by 1 and resets `currentDayOfYear` to 1.

## Design Concepts

**Stochasticity:** The model incorporates stochasticity in several processes. For example, the `set-parameters` procedure uses random distributions to assign values to the model parameters when the experiment type is set to "random". The `generate-annual-precipitation` procedure uses random-normal distributions to determine the parameters of the double logistic curve used to generate the annual precipitation pattern.

**Emergence:** The daily weather patterns emerge from the interaction of the various parametric equations and stochastic processes in the model. The seasonal patterns of temperature, precipitation, CO2 concentration, and solar radiation are not explicitly defined but emerge from the underlying equations.

**Sensing:** There is not a sensing component in this model.

**Adaptation:** There is not an adaptation component in this model.

**Objectives:** There are no agent objectives in this model.

**Learning:** There is no agent learning in this model.

**Prediction:** There is not a prediction component in this model.

**Collectives:** There are no agent collectives in this model.

**Observation:** The model allows for the observation of daily weather variables, including temperature, precipitation, CO2 concentration, solar radiation, and reference evapotranspiration. These variables can be used as inputs for other simulation models.

**Initialisation:** The model initializes the state variables at the beginning of the simulation. The `setup` procedure sets the `currentDayOfYear` variable to 1, generates the initial weather conditions using the `update-weather` procedure, and resets the ticks counter.

**Input Data:** The model does not use external input data. Instead, it relies on parametric equations and stochastic processes to generate the daily weather variables. However, the model can be calibrated by adjusting the model parameters to match historical weather data.

## Initialisation

The `setup` procedure initializes the model. It calls `set-constants` to define constant values such as `yearLengthInDays`. It calls `set-parameters` to set the model parameters based on the experiment type (user-defined or random). It sets `currentDayOfYear` to 1 and calls `update-weather` to generate the initial weather conditions. Finally, it resets the ticks counter.

## Input data

The model relies on parametric equations and stochastic processes to generate the daily weather variables. The parameters for these equations are set either by the user or through random distributions.

Here is a full account of the constants and modified parameters:

**Constants** are fixed values that do not change during the simulation. They represent fundamental properties or conditions of the environment being modeled. In the Weather model, the constants include: \* `yearLengthInDays`: Specifies the number of days in a year, set to 365. \* `albedo`: Represents the canopy reflection or albedo of a hypothetical grass reference crop, with a default value of 0.23. Albedo is a crucial factor in determining the net solar radiation. \* `elevation`: Indicates the elevation above sea level in meters, which could potentially be converted into a patch variable for more spatially explicit modeling.

**Modified parameters** are variables that can be adjusted to represent different environmental conditions or scenarios. These parameters influence the behavior and outputs of the model and can be modified through the user interface or experimental setup. They can be further grouped based on the specific aspect of weather they govern:

-   **Temperature (°C)**:
    -   `temperature_annualMaxAt2m`: The annual maximum daily average temperature at 2 meters.  
    -   `temperature_annualMinAt2m`: The annual minimum daily average temperature at 2 meters.  
    -   `temperature_meanDailyFluctuation`: The average daily fluctuation in temperature.  
    -   `temperature_dailyLowerDeviation`: Deviation of daily temperature to the lower side.  
    -   `temperature_dailyUpperDeviation`: Deviation of daily temperature to the upper side.  
-   **Precipitation (mm)**:
    -   `precipitation_yearlyMean`: The mean annual precipitation.  
    -   `precipitation_yearlySd`: The standard deviation of annual precipitation.  
    -   `precipitation_dailyCum_plateauValue_yearlyMean`: The mean plateau value of the annual cumulative precipitation curve.  
    -   `precipitation_dailyCum_plateauValue_yearlySd`: The standard deviation of the plateau value of the annual cumulative precipitation curve.  
    -   `precipitation_dailyCum_inflection1_yearlyMean`: Mean day of year in which the first logistic curve has its maximum slope.  
    -   `precipitation_dailyCum_inflection1_yearlySd`: Standard deviation of the day of year in which the first logistic curve has its maximum slope.  
    -   `precipitation_dailyCum_inflection2_yearlyMean`: Mean day of year in which the second logistic curve has its maximum slope.  
    -   `precipitation_dailyCum_inflection2_yearlySd`: Standard deviation of the day of year in which the second logistic curve has its maximum slope.  
    -   `precipitation_dailyCum_rate1_yearlyMean`: Mean maximum rate or slope increase of the first logistic curve.  
    -   `precipitation_dailyCum_rate1_yearlySd`: Standard deviation of the maximum rate or slope increase of the first logistic curve.  
    -   `precipitation_dailyCum_rate2_yearlyMean`: Mean maximum rate or slope increase of the second logistic curve.  
    -   `precipitation_dailyCum_rate2_yearlySd`: Standard deviation of the maximum rate or slope increase of the second logistic curve.  
    -   `precipitation_dailyCum_nSamples_yearlyMean`: Mean number of random samples or steps created in the cumulative curve.  
    -   `precipitation_dailyCum_nSamples_yearlySd`: Standard deviation of the number of random samples in the cumulative curve.  
    -   `precipitation_dailyCum_maxSampleSize_yearlyMean`: Mean maximum length of samples or step plateaus in the cumulative curve.  
    -   `precipitation_dailyCum_maxSampleSize_yearlySd`: Standard deviation of the maximum length of samples or step plateaus in the cumulative curve.  
-   **CO2 (ppm)**:
    -   `CO2_annualMin`: The annual minimum CO2 concentration.  
    -   `CO2_annualMax`: The annual maximum CO2 concentration.  
    -   `CO2_meanDailyFluctuation`: The average daily fluctuation in CO2 concentration.  
-   **Solar radiation (MJ/m2)**:
    -   `solar_annualMax`: The annual maximum solar radiation.  
    -   `solar_annualMin`: The annual minimum solar radiation.  
    -   `solar_meanDailyFluctuation`: The average daily fluctuation in solar radiation.  

## Procedures

**`gen-temperature`**: This procedure calculates the average temperature for the current day, using the `gen-annual-sinusoid-with-fluctuation` procedure.

**`gen-CO2`**: This procedure calculates the average CO2 concentration for the current day, using the `gen-annual-sinusoid-with-fluctuation` procedure.

**`gen-solar-radiation`**: This procedure calculates the solar radiation for the current day, using the `gen-annual-sinusoid-with-fluctuation` procedure.

**`update-precipitation`**: This procedure determines the precipitation for the current day by accessing a pre-generated annual precipitation series. The annual precipitation series is generated at the beginning of each model run, using a double logistic curve whose parameters are randomized based on yearly means and standard deviations.

**`get-ETr`**: This procedure calculates the reference evapotranspiration (`ETr`) based on the Penman-Monteith equation. It uses the daily temperature, solar radiation, and other parameters to estimate the `ETr` value.

**`gen-annual-sinusoid-with-fluctuation`:** Generates a temperature value for the current day, using a sine wave function to simulate the annual temperature cycle, with added daily fluctuation.

**`gen-annual-sinusoid`**: Generates a value based on a sine wave function, simulating an annual cycle.

**`gen-cumulative-precipitation-of-year`**: Simulates the cumulative proportion of yearly precipitation using a double logistic curve.

**`gen-annual-double-logistic-curve`**: Generates a double logistic curve, used as a proxy for the year series of daily cumulative precipitation.

**`get-increments-from-cumulative-curve`**: Derives daily proportions of yearly precipitation from the cumulative precipitation curve.

**`get-incremets-from-curve`**: Calculates the increments or derivatives corresponding to a given curve.

**`clampMinMax`**: Constrains a value within a specified minimum and maximum range.
