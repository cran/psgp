#include "astonGeostats.h"
#include <iostream>

/**
 * Function to convert metadata indices (1 to N) to C indices (0 to N-1)
 */
int minusOne(int n) { return n-1; }


void learnParameters(int numObs, double *xData, double *yData,  
        int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, 
        double* psgpParameters)
{
	vec rndNums = itpp::randu(numObs);
	ivec rndIdx = itpp::sort_index(rndNums);
	
	// Limit to MAX_OBS observation for parameter estimation
	// so that computations remain timely
	int nObsMax = min(numObs, MAX_OBS);

    // Set number of active points to max(MAX_ACTIVE_POINTS, size(obs))
    int n_active = min(MAX_ACTIVE_POINTS, nObsMax);

	
	mat X = itpp::zeros(nObsMax, 2);
	vec y = itpp::zeros(nObsMax);
	Vec<LikelihoodType *> noiseModels;

	// need to make this use a factory method from IT++
	// this would make things much more efficient
	for(int i=0; i < nObsMax; i++) {
		X(i, 0) = xData[rndIdx[i]];
		X(i, 1) = xData[rndIdx[i] + numObs];
		y(i) = yData[rndIdx[i]];
	}

	// RB: get variogram model parameters as a starting point
	// variogramParameters[0] is the model ID - we can ignore it.
	double range  = psgpParameters[1];
	double sill   = psgpParameters[2];
	double nugget = psgpParameters[3];
	double bias   = abs(1.0 / mean(y));    // Empirical estimation of bias
	
	// Make sure everything is fine - the variogram estimation can  
	// give invalid parameters (negative range...) - if not, revert
	// to some first guess from the data.
	if ( range <= 0.0 || sill <= 0.0 || nugget <= 0.0)
	{
	    cout << "Invalid variogram parameters: either the range, sill or nugget" << endl;
	    cout << "is negative or zero. Reverting to defaults." << endl;
	    
	    double r1=abs(min(max(X, 1) - min(X, 1)));
        double r2=abs(min(max(X, 2) - min(X, 2)));
        
        range = 0.25 * ((r1+r2) / 2.0);
        sill = abs(variance(y));
        nugget = 0.5 * sill;
	}
	
	// Covariance function
	ExponentialCF kernel1(range, sill);
	Matern5CF     kernel2(range, sill);
	ConstantCF    kernel3(bias);
	
	// Kernel made of several covariance functions
	SumCovarianceFunction covfKernel(kernel1);
	covfKernel.addCovarianceFunction(kernel2);
	covfKernel.addCovarianceFunction(kernel3);

	// Final covariance function is kernel + bias + white noise
	WhiteNoiseCF  covfNugget(nugget);
	SumCovarianceFunction covSum(covfKernel);
	covSum.addCovarianceFunction(covfNugget);
	
	cout << endl << "===Starting Parameters===========================" <<endl;
	covSum.displayCovarianceParameters();

	// By default, we use PSGP to estimate parameters, but there is
	// the option to use a GP instead.
	if(PARAMETER_ESTIMATION_USING_GP) 
	{
	    // Use a GP to learn parameters
	    GaussianProcess gp(2, 1, X, y, covSum);
	
    	SCGModelTrainer gpTrainer(gp);
    
    	gpTrainer.setAnalyticGradients(true);
    	gpTrainer.setCheckGradient(false);
    	gpTrainer.Train(50);
	}
	else
	{
	    int numModels = numMetadata;       // Total number of likelihood models
	    GaussianLikelihood* defaultLikelihood;
	    Vec<LikelihoodType*> likelihoodModels(numModels);        // Vector of likelihood models
	    ivec iLikelihoodModel(nObsMax);
	    
	    // PARSE METADATA IF AVAILABLE, OTHERWISE USE DEFAULT LIKELIHOOD MODEL
	    if(numMetadata == 0)
	    {
	        cout << "No noise model specified" << endl;
	        cout << "Defaulting to GAUSSIAN with variance " << (nugget * LIKELIHOOD_NUGGET_RATIO) << endl;
	        defaultLikelihood = new GaussianLikelihood(nugget * LIKELIHOOD_NUGGET_RATIO);
	    }
	    else
	    {
	        cout << "Observation error characteristics specified. " << endl;
	        cout << "Building error models from metadata table." << endl;

	        Vec<string> modelName(numModels);      // Name of likelihood model for each obs
	        Vec<vec>    modelParams(numModels);    // Params of likelihood model for each obs
	        
	        // Shift indexes (starting from 1 in data, from 0 in ITPP)
	        ivec iLikelihoodModelAllObs = apply_function(minusOne, ivec(errPtr, numObs));
	        
	        // Retain only the indexes for the observations considered
	        iLikelihoodModel = iLikelihoodModelAllObs(rndIdx(0,nObsMax-1));

	        // Extract the noise models and their parameters for all observations
	        parseMetadata(metaDataTable, modelName, modelParams);

	        // Build vector of likelihod models 
	        buildLikelihoodVector(modelName, modelParams, likelihoodModels, 
	                              nugget * LIKELIHOOD_NUGGET_RATIO);
	    }
	    
	    // LEARN PARAMETERS USING PSGP AND LIKELIHOOD MODEL(S) GENERATED ABOVE
	    PSGP psgp(X, y, covSum, n_active, NUM_SWEEPS_CHANGING, NUM_SWEEPS_FIXED);
	            
	    if (numMetadata == 0) {
	        psgp.computePosterior(*defaultLikelihood);
	    }
	    else {
	        psgp.computePosterior(iLikelihoodModel, likelihoodModels);
	    }
	    
	    
	    // Use approximate objective function (more stable)
	    psgp.setLikelihoodType(Approximate);
	    
	    SCGModelTrainer gpTrainer(psgp);           // SCG optimisation method
	    gpTrainer.setAnalyticGradients(true);
	    gpTrainer.setCheckGradient(false);
	    
	    cout << endl << "Finding optimal parameters" << endl;
	    for (int i=0; i<PSGP_PARAM_ITERATIONS; i++)
	    {
	        gpTrainer.Train(PSGP_SCG_ITERATIONS);
	        psgp.recomputePosterior();
	    }
	}
	
	cout << endl << "===Final Parameters===========================" <<endl;
	covSum.displayCovarianceParameters();

	// Copy final parameters over to psgpParameters
	vec finalParams = covSum.getParameters();
	for(int i=0; i<finalParams.length(); i++)
	{
	    *psgpParameters++ = finalParams(i);
	}
		
	// Add padding zeros (remember psgpParameters has fixed size and is
	// likely to be bigger than we need)
	for(int i=finalParams.length(); i<NUM_PSGP_PARAMETERS; i++)
	{
	    *psgpParameters++ = 0.0;
	}
	
}


/**
 * Prediction using Projected Sequential Gaussian Processes at a set of specified locations
 * 
 * int numObs               Number of observations
 * int numPred              Number of predictions
 * double *xData            Array of observed locations
 * double *yData            Array of observed values
 * double *xPred            Array of prediction locations
 * int numMetadata          Size of metadata table (i.e. number of likelihood models)
 * int *errPtr              Array of likelihood model indexes (e.g. 1 per observation)
 * int *sensorPtr           Array of sensor indexes (e.g. 1 per observation)
 * char **metaDataTable,    Metadata array of strings, each string containing the description of a 
 *                          likelihood model in the form "<DISTRIBUTION_NAME>,<PARAM_1>,...,<PARAM_N>" 
 *
 * double *meanPred,        Array of prediction mean values                
 * double *varPred,         Array of prediction variances
 * 
 * double *range,           PSGP range parameter                
 * double *sill,            PSGP sill parameter
 * double *nugget,          PSGP nugget parameter
 * double *bias,            PSGP bias parameter 
 * int model                ? (Unused)
 */
void makePredictions(int numObs, int numPred, double *xData, double *yData,  double *xPred, 
                     int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, 
                     double *meanPred, double *varPred, double* psgpParameters)
{
    LikelihoodType *defaultLikelihood;
    Vec<LikelihoodType*> likelihoodModels(numObs);        // Vector of likelihood models
	ivec iLikelihoodModel(numObs);                        // Index of likelihood model for each obs
	
	// Convert arguments to IT++ vectors/matrices
	vec rndNums = itpp::randu(numObs);
	ivec rndIdx = itpp::sort_index(rndNums);

	mat X = reshape(vec(xData,2*numObs),numObs,2);
	vec y = vec(yData,numObs);
	
	mat Xpred = reshape(vec(xPred,2*numPred),numPred,2);
	vec predY = itpp::zeros(numPred);
	vec predVar = itpp::zeros(numPred);
	
	//-------------------------------------------------------------------------
    // COVARIANCE FUNCTION
	// 
    // The covariance function parts are created with dummy parameter values.
	// The full set of parameters is then overwritten by properly estimated    
	// values.
    //-------------------------------------------------------------------------
	double dummy_range  = 1.0;   // This is a dummy value - overwritten below.
	double dummy_sill   = 1.0;   // This is a dummy value - overwritten below.
	double dummy_bias   = 1e-3;  // This is a dummy value - overwritten below.
	double dummy_nugget = 1e-2;  // This is a dummy value - overwritten below.
	
	// Covariance function
    ExponentialCF kernel1(dummy_range, dummy_sill);
    Matern5CF     kernel3(dummy_range, dummy_sill);
    ConstantCF    kernel5(dummy_bias);
    
    // Kernel made of several covariance functions
    // This is the kernel used for prediction
    SumCovarianceFunction covfKernel(kernel1);
    covfKernel.addCovarianceFunction(kernel3);
    covfKernel.addCovarianceFunction(kernel5);

    // Final covariance function is kernel + white noise
    WhiteNoiseCF  covfNugget(dummy_nugget);
    SumCovarianceFunction covSum(covfKernel);
    covSum.addCovarianceFunction(covfNugget);
	
    // Overwrite current dummy parameters with correct, estimated ones
    int nparams = covSum.getNumberParameters();
    covSum.setParameters( vec(psgpParameters, nparams) );
	
    // Nugget is the last parameter
    double nugget = covSum.getParameter(nparams-1);
    
	
	//-------------------------------------------------------------------------
	//
	// PARSE METADATA IF AVAILABLE, OTHERWISE USE DEFAULT LIKELIHOOD MODEL
	//
	//-------------------------------------------------------------------------
	if(numMetadata == 0)
	{
		cout << "No noise model specified" << endl;
		cout << "Defaulting to GAUSSIAN with variance = " << (nugget * LIKELIHOOD_NUGGET_RATIO) << endl;
		defaultLikelihood = new GaussianLikelihood(nugget * LIKELIHOOD_NUGGET_RATIO);
	}
	else
	{
	    cout << "Observation error characteristics specified. " << endl;
	    cout << "Building error models from metadata table." << endl;
	 
	    int numModels = numMetadata;              
	    Vec<string>  modelName(numModels);     // Names of likelihood models
	    Vec<vec>     modelParams(numModels);   // Params of likelihood models
	    
	    // Shift indexes (starting from 1 in data, from 0 in ITPP)
	    iLikelihoodModel = apply_function(minusOne, ivec(errPtr, numObs));
	    
	    // Make sure we have same number of locations and error models
	    // assert(iLikelihoodModel.length() == numMetadata);
	    // assert(iLikelihoodModel.length() == X.rows());
	    
	    // Extract the noise models and their parameters for all observations
	    parseMetadata(metaDataTable, modelName, modelParams);
	    
	    // Build vector of likelihod models 
	    buildLikelihoodVector(modelName, modelParams, likelihoodModels, 
	                          nugget * LIKELIHOOD_NUGGET_RATIO);
	}
	
	
	//-------------------------------------------------------------------------
	//
	// INITIALISE PSGP
	//
	//-------------------------------------------------------------------------
	
	// Set number of active points to max(MAX_ACTIVE_POINTS, size(obs))
	int n_active = min(MAX_ACTIVE_POINTS, numObs);
	
	// SequentialGP ssgp(2,1,n_active,X,y,covSum,1);
	PSGP ssgp(X, y, covSum, n_active, NUM_SWEEPS_CHANGING, NUM_SWEEPS_FIXED);

	cout << "Computing posterior..." << endl;
	if (numMetadata == 0) {
	    ssgp.computePosterior(*defaultLikelihood);
	}
	else {
	    ssgp.computePosterior(iLikelihoodModel, likelihoodModels);
	}
	
	if ( !USING_CHUNK_PREDICTION )
	{
	    cout << "Predicting..." << endl;
	    ssgp.makePredictions(predY, predVar, Xpred, covfKernel);
	}
	else
	{
	    int startVal = 0;
	    int chunkSize = CHUNK_SIZE;
	    int endVal = chunkSize - 1;

	    if(endVal > numPred)
	    {
	        endVal = numPred - 1;
	        chunkSize = endVal - startVal + 1;
	    }

	    while(startVal < numPred)
	    {
	        cout << "  Predicting chunk [" << startVal << ":" << endVal << "/" << numPred << "]" << endl;
	        
	        mat XpredChunk = Xpred.get_rows(startVal, endVal);

	        // Predicted mean and variance for data chunk
	        vec predYChunk(chunkSize);
	        vec predVarChunk(chunkSize);

	        ssgp.makePredictions(predYChunk, predVarChunk, XpredChunk, covfKernel);

	        predY.replace_mid(startVal, predYChunk);
	        predVar.replace_mid(startVal, predVarChunk);
	        
	        startVal = endVal + 1;
	        endVal = endVal + chunkSize;
	        
	        if(endVal >= numPred)
	        {
	            endVal = numPred - 1;
	            chunkSize = endVal - startVal + 1;
	        }
	    }
	}

	// should use an IT++ factory for vectors
	for(int i=0; i < numPred; i++)
	{
		meanPred[i] = predY(i);
		varPred[i]  = predVar(i);
	}

	cout << "PSGP used the following parameters:" << endl;
	covSum.displayCovarianceParameters();

	cout << "Done." << endl;
}


/**
 * Parse the metadata table.
 * 
 * PARAMETERS:
 * 
 *  metadataTable:  a array of C strings, each element being a textual
 *                  representation of a distribution (noise model) and
 *                  its parameters (e.g. "GAUSSIAN,0.0,1.0" for the normal
 *                  distribution)
 *  
 *  modelIndex:     a vector of indexes giving, for each observation, 
 *                  the index of the row in the metadata table that
 *                  contains the corresponding noise model data.
 * 
 * RETURNS:
 * 
 *  vecModelName:   a std::vector of distribution names (strings). These
 *                  correspond to the first part of the metadata (e.g.
 *                  "GAUSSIAN")
 * 
 *  vecModelParams: a std::vector of distribution parameters (strings).
 *                  These correspond to the parameters (as text) of the 
 *                  distribution (e.g. "0.0,1.0")
 * 
 *  modelIndex:     Observations without a valid noise model have their index
 *                  set to -1 in this vector, for further processing in the
 *                  calling function. 
 */  
void parseMetadata(char** metadataTable, Vec<string> &modelNames, Vec<vec> &modelParams)
{
    // For each observation
    for(int iModel = 0; iModel < modelNames.length(); iModel++)
    {
        // modelSet(iObs) = false;
        
        // Get index of corresponding noise model
        // int iNoiseModel = modelIndex(iObs);
        
        // Retrieve corresponding string in metadata table
        // Should be something like "GAUSSIAN,0.0,0.123", i.e. a list of string
        // tokens where the first is the distribution identifier and the following
        // are the values for the parameters (number of these depends on the distribution
        // type, e.g. Gaussian has 2 parameters: mean and variance)
        // cout << "metadataTable[" << iNoiseModel << "]: " << metadataTable[iNoiseModel] << endl;
        string strNoiseModel(metadataTable[iModel]);
        
        // Split string in 2 w.r.t. to first (non-trailing) comma delimiter
        std::vector<string> tokens;
        itppext::tokenise(strNoiseModel, tokens, ", ", 1);

        // Retrieve distribution name
        modelNames[iModel] = tokens[0];
        
        // Use IT++ vector initialisation from string - if this does not work, 
        // i.e. string cannot be converted to double, flag the noise model as
        // invalid by setting the modelIndex to -1
        try {
            modelParams[iModel] = vec(tokens[1]);
        }
        catch (std::runtime_error e) {
            cout << "** Error in metadata parsing for noise model " << iModel << ":" << endl;
            cout << "   Invalid parameter string \"" << tokens[1] << "\"" << endl;
            cout << "   Parameter string must be a sequence of numeric values" << endl;
            cout << "   separated by commas, e.g. \"1.23,4,5.6,78.9\"" << endl;
         
            // Flag observation as not having a valid model
            modelNames[iModel] = INVALID_MODEL_NAME;
        }
    }
    
}

/**
 * Return a vector of likelihood models based on array of model names
 * and array of model parameters.
 * Observations with an invalid noise model are given a default Gaussian 
 * likelihood with variance being:
 * - the average of the other Gaussian model's variances,
 *   if some valid Gaussian models were found
 * - a default value (defaulVar) if no valid Gaussian model was found 
 */
void buildLikelihoodVector(Vec<string> modelName, Vec<vec> modelParams, 
                           Vec<LikelihoodType*> &likelihoodModels, double defaultVar)
{
    int numModels  = modelName.length();
    
    int countGaussianModels = 0;   // The number of Gaussian noise models correctly parsed
    double avgErrMean = 0.0, avgErrVar = 0.0;  // The average mean and variance of these models

    for(int iModel=0; iModel<numModels; iModel++) 
    {
        // Get parameter vector for current model
        vec params = modelParams[iModel];

        // If a valid model name/params exist for this obs, 
        // create likelihood model
        if (modelName[iModel] == "GAUSSIAN") 
        {
            // TODO: At the moment, the GaussianLikelihood takes only the variance parameter
            // Change to GaussianLikelihood(mean, variance)
            // Also, add support for LikelihoodType(vec params) to make it more generic
            likelihoodModels[iModel] = new GaussianLikelihood( params(1) );
            
            // TODO: uncomment line below once the mean is taken into account in GaussianLikelihood
            // avgMean += distParams(0);
            avgErrVar += modelParams[iModel](1);
            countGaussianModels++;
        }
        /* Use template below to add other noise models 
                ---
                else if (modelName[iObs] == "MY_MODEL")
                {
                    myLikelihoods.push_back( MyLikelihoodType( params(0), ..., params(n) ) );
                    likelihoodModels[iModel] = new MyLikelihoodType( params(0), ..., params(n) );
                }
                ---
         */
        else
        {
            cout << "Unrecognized observation noise model: " << modelName[iModel] << 
            "for observation " << iModel << endl;

            // Flag observation as having an invalid noise model 
            modelName[iModel] = INVALID_MODEL_NAME;
        }
    }

    // If we found some valid Gaussian noise models, infer parameters for
    // default likelihood from theirs (take average mean and average variance).
    // Should hopefuly be slightly better than our first guess.
    if (countGaussianModels > 1) {
        avgErrMean = avgErrMean/countGaussianModels;
        avgErrVar  = avgErrVar/countGaussianModels;
    }
    // If no valid Gaussian model is found, revert to the default Gaussian 
    // likelihood
    else
    {
        avgErrVar = defaultVar;
    }

    cout << "Observations without a valid noise model are given a default" << endl;
    cout << "  Gaussian noise model with (mean, variance) = (" << avgErrMean; 
    cout << "," << avgErrVar << ")" << endl;

    // Attribute all observations with an invalid noise model the default 
    // (Gaussian) likelihood model
    for(int iModel=0; iModel<numModels; iModel++) 
    {
        if(modelName[iModel] == INVALID_MODEL_NAME) 
        {
            // likelihoodModels[iLikelihoodModel(iModel)] = defaultLikelihood;
            likelihoodModels[iModel] = new GaussianLikelihood(avgErrVar);
        }
    }
}

