estimateParameters.psgp = function(object,...) {

	origObs = object$observations
	
	rotated = FALSE
	if (object$params$doAnisotropy) {
  		object = estimateAnisotropy(object) 
	  	#rotate Data
  	if (object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1"))
			object$observations=rotateAnisotropicData(object$observations,object$anisPar)
  	rotated = TRUE
	}  

	#if (require(astonGeostats)) 
	object = learnParameters(object)
	if (rotated) 
		object$observations = origObs
	object
}


spatialPredict.psgp = function(object,...) {

dots = list(...)
variogramParameters = object$variogramModel

vario=array()

# variogram type
# Gau - 1
# Exp - 2
#

vario[1] = as.integer(object$variogramModel$model[2])
#if(object$variogramModel$model[2] == "Gau") vario[1]=1
#if(object$variogramModel$model[2] == "Exp") vario[1]=2
vario[2]=object$variogramModel$range[2]
vario[3]=object$variogramModel$psill[2]
vario[4]=object$variogramModel$psill[1]
vario[5]=object$variogramModel$beta[1]

rotated = FALSE
if (object$params$doAnisotropy && object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1")){
  objTemp = object
  object$observations = rotateAnisotropicData(object$observations, object$anisPar)
  object$predictionLocations = rotateAnisotropicData(object$predictionLocations, object$anisPar)
  rotated = TRUE
}

#if (require(astonGeostats)) {
  p = makePrediction(object, vario)
  object$predictions = SpatialPointsDataFrame(object$predictionLocations,
    data = data.frame(var1.pred = unlist(p[1]),var1.var=unlist(p[2])))
  nsim = ifelse("nsim" %in% names(dots),dots$nsim,0) 
  if (nsim > 0) {
    nmax = object$params$nmax
    object$predictions@data = cbind(object$predictions@data,krige(object$formulaString,object$observations, 
           object$predictionLocations,object$variogramModel,nsim=nsim,nmax = nmax,debug.level = object$params$debug.level)@data)
  }
  if (rotated) {
    object$observations = objTemp$observations
    object$predictionLocations = objTemp$predictionLocations
    object$predictions@coords = coordinates(object$predictionLocations)
    object$predictions@bbox = bbox(object$predictionLocations)
    proj4string(object$predictions) = proj4string(object$predictionLocations)
  }
  names(object$predictions)[1:2] = c("var1.pred","var1.var")
#}
object
}


summary.psgp = function(object, ...) {
  summaryIntamap(object, ...)
}



plot.psgp = function(x, ...) {
  plotIntamap(x,  ...)
}
