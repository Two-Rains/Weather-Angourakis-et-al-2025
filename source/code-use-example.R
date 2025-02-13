
source("source/weatherModel.R")

SEED <- 0
YEAR_LENGTH <- 365 # ignoring leap year adjustment
NUM_YEARS <- 5
NUM_DAYS <- NUM_YEARS * YEAR_LENGTH

weather_model <- initialise_weather_model(seed = SEED, year_length = YEAR_LENGTH)

weather_model <- run_weather_model(weather_model, num_years = NUM_YEARS)
