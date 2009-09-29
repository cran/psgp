makePrediction <- function(object, vario)
{

inputs = object$observations
pred = object$predictionLocations

# put data into an easy parseable format for the backend C++ code
x = coordinates(inputs)
y = as.vector(inputs$value)

### error variance vector
e = as.vector(inputs$var)

tx = coordinates(pred)

### error variance vector
e = as.vector(inputs$value)

# put obsError and sensorID into vectors
observationError = as.integer(inputs$oeid)
sensorModel = as.integer(inputs$sensor)

# vector of strings of meta data
metaData = object$obsChar


r <- .Call("predict", x, y, e, tx, vario, observationError, sensorModel, metaData,
	PACKAGE = "psgp")

}

