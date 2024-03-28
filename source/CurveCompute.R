# load datasets ################################################################################

# inputFiles <- paste0("inputData/", list.files(path = "inputData"))
# 
# weatherData <- data.frame()
# 
# for (i in 1:length(inputFiles))
# {
#   tempData <- read.csv(inputFiles[i], skip = 18)
#   tempData <- tempData[tempData$YEAR %in% 1984:2007,]
#   
#   weatherData <- rbind(weatherData, tempData)
# }
# 
# rm(tempData)
# 
# # use Latitude to distinguish sites and create a new variable using site names:
# weatherData$Site <- rep(NA, nrow(weatherData))
# for (i in 1:nrow(weatherData))
# {
#   if (floor(weatherData$LAT[i]) == 29) { weatherData$Site[i] <- "Rakhigarhi" }
#   if (floor(weatherData$LON[i]) == 104) { weatherData$Site[i] <- "Irkutsk" }
#   if (floor(weatherData$LAT[i]) == -43) { weatherData$Site[i] <- "Hobart" }
#   if (floor(weatherData$LAT[i]) == 21) { weatherData$Site[i] <- "Pearl Harbor" }
#   if (floor(weatherData$LAT[i]) == -24) { weatherData$Site[i] <- "Sao Paulo" }
#   if (floor(weatherData$LON[i]) == 0) { weatherData$Site[i] <- "Cambridge" }
#   if (floor(weatherData$LAT[i]) == -23) { weatherData$Site[i] <- "Windhoek" }
# }

weatherData <- read.csv("weatherData.csv")

# compute MSE between simulated and one data curve ################################################################################

curveCompute <- function(
    locationName,
    year,
    plateauValue, 
    inflection1, 
    rate1, 
    inflection2, 
    rate2, 
    nSamples, 
    maxSampleSize, 
    seed = 0
)
{
  
    dataPrecipitationYear <- weatherData$PRECTOT[weatherData$Site == locationName & 
        weatherData$YEAR == year]
    
    simulatedPrecipitationYear <- getPrecipitationOfYear(
      plateauValue = plateauValue, 
      inflection1 = inflection1, 
      rate1 = rate1, 
      inflection2 = inflection2, 
      rate2 = rate2, 
      yearLengthInDays = length(dataPrecipitationYear), 
      nSamples = nSamples, 
      maxSampleSize = maxSampleSize, 
      annualSum = sum(dataPrecipitationYear),
      seed = seed
    )
    
    targetCurve <- getCumulativePrecipitationOfYear(dataPrecipitationYear)
    simulatedCurve <- getCumulativePrecipitationOfYear(simulatedPrecipitationYear)
    
    MSE <- mean((targetCurve - simulatedCurve)^2)
    
    return(MSE)
}

# compute MSE between simulated and the average of several data curves ################################################################################

curveComputeMultipleYears <- function(
    locationName,
    plateauValue, 
    inflection1, 
    rate1, 
    inflection2, 
    rate2, 
    nSamples, 
    maxSampleSize, 
    seed = 0
)
{
  
  dataPrecipitationYear <- weatherData$PRECTOT[weatherData$Site == locationName]
  
  simulatedPrecipitationYear <- getPrecipitationOfYear(
    plateauValue = plateauValue, 
    inflection1 = inflection1, 
    rate1 = rate1, 
    inflection2 = inflection2, 
    rate2 = rate2, 
    yearLengthInDays = length(dataPrecipitationYear), 
    nSamples = nSamples, 
    maxSampleSize = maxSampleSize, 
    annualSum = sum(dataPrecipitationYear),
    seed = seed
  )
  
  # calculate all curves per year
  listOfYearCurves <- list()
  for (year in levels(factor(weatherData$YEAR)))
  {
    tempCurve <- getCumulativePrecipitationOfYear(dataPrecipitationYear[weatherData$YEAR == year])
    
    listOfYearCurves <- append(listOfYearCurves, tempCurve)
  }
  # find average curve
  targetCurve <- Reduce("+",listOfYearCurves)/length(listOfYearCurves)
  
  simulatedCurve <- getCumulativePrecipitationOfYear(simulatedPrecipitationYear)
  
  MSE <- mean((targetCurve - simulatedCurve)^2)
  
  return(MSE)
}

# run ################################################################################

# delta <-
#   curveCompute(
#     locationName = locationName,
#     year = year,
#     plateauValue = plateauValue,
#     inflection1 = inflection1,
#     rate1 = rate1,
#     inflection2 = inflection2,
#     rate2 = rate2,
#     nSamples = nSamples,
#     maxSampleSize = maxSampleSize,
#     seed = seed
#   )

delta <-
  curveComputeMultipleYears(
    locationName = locationName,
    plateauValue = plateauValue,
    inflection1 = inflection1,
    rate1 = rate1,
    inflection2 = inflection2,
    rate2 = rate2,
    nSamples = nSamples,
    maxSampleSize = maxSampleSize,
    seed = seed
  )
