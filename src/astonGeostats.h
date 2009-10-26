#ifndef ASTONGEOSTATS_H
#define ASTONGEOSTATS_H


#include "R.h"
#include "Rmath.h"

#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "gaussianProcesses/ForwardModel.h"
#include "gaussianProcesses/PSGP.h"

#include "gaussianProcesses/GaussianProcess.h"
#include "optimisation/Optimisable.h"
#include "optimisation/SCGModelTrainer.h"
#include "optimisation/QuasiNewtonModelTrainer.h"
#include "optimisation/CGModelTrainer.h"
#include "optimisation/GDModelTrainer.h"
#include "optimisation/ModelTrainer.h"

#include "likelihoodModels/LikelihoodType.h"
#include "likelihoodModels/ExponentialSampLikelihood.h"
#include "likelihoodModels/GaussianLikelihood.h"

#include "covarianceFunctions/CovarianceFunction.h"
#include "covarianceFunctions/GaussianCF.h"
#include "covarianceFunctions/ExponentialCF.h"
#include "covarianceFunctions/WhiteNoiseCF.h"
#include "covarianceFunctions/Matern3CF.h"
#include "covarianceFunctions/Matern5CF.h"
#include "covarianceFunctions/NeuralNetCF.h"
#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/LogTransform.h"
#include "covarianceFunctions/IdentityTransform.h"
#include "covarianceFunctions/NegLogSigmoidTransform.h"
#include "covarianceFunctions/ConstantCF.h"

#include "itppext/itppext.h"
#include "io/csvstream.h"

#include <cassert>
#include <string>
#include <algorithm>
#include <stdexcept>

//-----------------------------------------------------------------------------
// CONSTANTS AND OTHER GENERAL PARAMETERS
//-----------------------------------------------------------------------------

// Max number of parameters for PSGP (this limit is set by the R code)
#define NUM_PSGP_PARAMETERS 16

// Maximum number of observations kept for param estimation
#define MAX_OBS 1000

// Maximum number of active points
#define MAX_ACTIVE_POINTS 400

// Likelihood to nugget ratio
#define LIKELIHOOD_NUGGET_RATIO 0.01

// Number of sweeps through data with changing/fixed active set
#define NUM_SWEEPS_CHANGING 1
#define NUM_SWEEPS_FIXED 1

// ID for invalid noise model (when using observation noise) 
#define INVALID_MODEL_NAME "INVALID_MODEL"

// Whether to use a GP instead of PSGP for parameter estimation
#define PARAMETER_ESTIMATION_USING_GP false

// Outer loops in parameter estimation for PSGP
#define PSGP_PARAM_ITERATIONS 3

// Inner loop (i.e. SCG iterations in each outer loop) for PSGP
#define PSGP_SCG_ITERATIONS 5

// Define whether to use full prediction (all data at once) or 
// split prediction (by chunks of data)
#define USING_CHUNK_PREDICTION true

// Size of prediction chunk (in number of observations)
#define CHUNK_SIZE 1000



using namespace std;
using namespace itpp;

// RB: These following been replaced by psgpInitParams (model was not used so removed too)
//     double *range, double *sill, double *nugget, double *bias 
//     int model
void learnParameters(int numObs, double *xData, double *yData,  
                     int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, 
                     double* psgpParams);

// RB: These have been replaced by psgpParams (model was not used so removed too)
//     double *range, double *sill, double *nugget, double *bias
//     int model
void makePredictions(int numObs, int numPred, double *xData, double *yData, 
                     double *xPred, int numMetadata, int *errPtr, int *sensorPtr, 
                     char **metaDataTable, double *meanPred, double *varPred, double* psgpParams);

void parseMetadata(char** metadataTable, Vec<string> &modelNames, Vec<vec> &modelParams);

void buildLikelihoodVector(Vec<string> modelName, Vec<vec> modelParams,  
                           Vec<LikelihoodType*> &likelihoodModels, double defaultVar);

int minusOne(int n);


#endif // ASTONGEOSTATS_H
