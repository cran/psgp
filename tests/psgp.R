rm(list = ls(all.names = TRUE))

# Suppress warnings and OpenMP messages
options(warn = -1)
options(save = "no")
options(digits=8)
Sys.unsetenv("KMP_DEVICE_THREAD_LIMIT")
Sys.unsetenv("KMP_ALL_THREADS")
Sys.unsetenv("KMP_TEAMS_THREAD_LIMIT")
Sys.unsetenv("OMP_THREAD_LIMIT")

library(automap)
library(psgp)

set.seed(13531)

data(meuse)
observations <- data.frame(x = meuse$x,y = meuse$y,value = log(meuse$zinc))
coordinates(observations) = ~x+y

predictionLocations <- spsample(observations, 50, "regular")

krigingObject <- createIntamapObject(
	observations = observations,
	predictionLocations = predictionLocations,
  formulaString = as.formula(value~1),
	params =  list(doAnisotropy = TRUE, thresh = quantile(observations$value,0.9)),
  outputWhat = list(mean=TRUE, variance=TRUE, excprob = 5.9, cumdistr = 5.9, 
		quantile = .1)
)
class(krigingObject) <- c("psgp")

checkSetup(krigingObject)
krigingObject <- preProcess(krigingObject)
krigingObject <- estimateParameters(krigingObject)
krigingObject <- spatialPredict(krigingObject)
krigingObject <- postProcess(krigingObject)

# Send predictions back to Java. 
summary(krigingObject$outputTable)
summary(krigingObject$observations)
summary(autoKrige(value~1,krigingObject$observations,predictionLocations)$krige_output)
autofitVariogram(value~1,krigingObject$observations)$var_model

# Restore original settings at the end
options(warn = 0)  # Restore warning level
