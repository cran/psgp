spatialPredict.psgp = function(object,...) {

    dots = list(...)
    
    # PSGP parameters are stored in the variogram model data frame
    # This is a hack to remain compatible with intamap. The first column
    # is discarded as it is for 
    # text and used to flag the values as incorrect (to ensure
    # this is not used as a variogram model). The other 9
    # columns contain the log parameters, which we put back into a vector
    params = object$variogramModel
    psgpLogParams = c(as.numeric(params[1,2:9]), as.numeric(params[2,2:9]))
    
    # vario=array()
    
    # variogram type
    # Gau - 1
    # Exp - 2
    #
    
    # vario[1] = as.integer(object$variogramModel$model[2])
    #if(object$variogramModel$model[2] == "Gau") vario[1]=1
    #if(object$variogramModel$model[2] == "Exp") vario[1]=2
    #vario[2]=object$variogramModel$range[2]
    #vario[3]=object$variogramModel$psill[2]
    #vario[4]=object$variogramModel$psill[1]
    #vario[5]=object$variogramModel$beta[1]
    
    
    rotated = FALSE
    if (object$params$doAnisotropy && object$anisPar$doRotation && all(as.character(object$formulaString[[3]])=="1")){
      objTemp = object
      object$observations = rotateAnisotropicData(object$observations, object$anisPar)
      object$predictionLocations = rotateAnisotropicData(object$predictionLocations, object$anisPar)
      rotated = TRUE
    }
    
    #if (require(astonGeostats)) {
      p = makePrediction(object, psgpLogParams)
      object$predictions = SpatialPointsDataFrame(object$predictionLocations,
        data = data.frame(var1.pred = unlist(p[1]),var1.var=unlist(p[2])))
      nsim = ifelse("nsim" %in% names(dots),dots$nsim,0) 
      if (nsim > 0) {
        nmax = object$params$nmax
        object$predictions = cbind(object$predictions,krige(object$formulaString,object$observations, 
               object$predictionLocations,object$variogramModel,nsim=nsim,nmax = nmax,debug.level = object$params$debug.level))
      }
      if (rotated) {
        object$observations = objTemp$observations
        object$predictionLocations = objTemp$predictionLocations
        object$predictions@coords = coordinates(object$predictionLocations)
        object$predictions@bbox = bbox(object$predictionLocations)
        proj4string(object$predictions) = proj4string(object$predictionLocations)
      }
      names(object$predictions) = c("var1.pred","var1.var")
    #}
    object
}
