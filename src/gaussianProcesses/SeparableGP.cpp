#include "SeparableGP.h"

/***
 * Creates a new instance of SeparableGP
 * 
 * @param inputs   A matrix of data inputs
 * @param outputs  A matrix of data outputs (evaluated at inputs)
 * @param cf_inputs The covariance function between inputs
 * @param cf_outputs The covariance function between outputs
 **/
SeparableGP::SeparableGP(mat _inputs, mat _outputs, CovarianceFunction& cf_inputs, CovarianceFunction& cf_outputs)
: ForwardModel(inputs.cols(), outputs.cols()), covFuncIn(cf_inputs), covFuncOut(cf_outputs)
{
    inputs = _inputs;
    outputs = _outputs;
}

/**
 * Destructor
 **/
SeparableGP::~SeparableGP()
{
}

/**
 * Make predictions at a set of new inputs. The covariance function between inputs can be specified. 
 * This is useful to make non-noisy predictions. The mean and variance predicted by the GP are returned. 
 * 
 * @param new_inputs The matrix of a new inputs
 * @param cf_inputs The covariance matrix between inputs 
 * @return gpmean The predicted mean
 * @return gpvar  The predicted variance  
 **/ 
void SeparableGP::makePredictions(vec &gpmean, mat &gpvar, mat new_inputs, CovarianceFunction &cf_inputs)
{
    mat Kf, Kx, k, ks;
    covFuncIn.computeSymmetric(Kf, inputs);                       // Kf = Sigma
    covFuncIn.computeSymmetric(ks, new_inputs);                   // k* = K(x*,x*)
    covFuncOut.computeSymmetric(Kx, outputs);                     // Kx = K(X,X)
    covFuncIn.computeCovariance(k, inputs, new_inputs);           // k  = K(X,x*)
    
    mat alpha = ls_solve_chol(Kx,outputs);                         // a = K(X,X)^{-1}*y
    mat gamma = backslash(chol(Kx), k);
    
    
    gpmean = cvectorize(k*alpha);
    gpvar  = kron(Kf,ks - gamma.transpose()*gamma);
}


/**
 * Make predictions at a set of new inputs. The current covariance function is used; if this covariance 
 * includes a noise term, the predictions will be noisy. The mean and variance predicted by the GP 
 * are returned. 
 * 
 * @param new_inputs The matrix of a new inputs
 * @return gpmean The predicted mean
 * @return gpvar  The predicted variance  
 **/
void SeparableGP::makePredictions(vec &gpmean, mat &gpvar, mat new_inputs)
{
   
}

//=============================================================================
//
// MODIFIERS
//
//=============================================================================

/**
 * Sets the covariance function between inputs
 *
 * @param cf A pointer to a valid covariance function object
 */
void SeparableGP::setCovarianceFunctionInputs(CovarianceFunction& cf)
{
    covFuncIn = cf;
}

/**
 * Sets the covariance function between outputs
 *
 * @param cf A pointer to a valid covariance function object
 */
void SeparableGP::setCovarianceFunctionOutputs(CovarianceFunction& cf) 
{
    covFuncOut = cf;
}

/**
 * Sets the parameter estimation method. 
 * 
 * @param method Can be MARGINAL_LOG_LIKELIHOOD or MARGINAL_LOG_PRIOR
 */
void SeparableGP::setParameterEstimationMethod(ParameterEstimationMethod method)
{
    paramEstimationMethod = method;
}

//=============================================================================
//
// Accessors
//
//=============================================================================

/**
 * Returns a pointer to the covariance function between inputs
 */
CovarianceFunction& SeparableGP::getCovarianceFunctionInputs()
{
    return covFuncIn;
}

/**
 * Returns a pointer to the covariance function between outputs
 */
CovarianceFunction& SeparableGP::getCovarianceFunctionOutputs()
{
    return covFuncOut;
}

/**
 * Returns the name of the current parameter estimation method
 */
string SeparableGP::getParameterEstimationMethodName() 
{
    return parameterEstimationMethodName[paramEstimationMethod];
}

//=============================================================================
//
// PARAMETER ESTIMATION METHODS (IMPLEMENTING OPTIMISABLE)
//
//=============================================================================
vec SeparableGP::getParametersVector() const
{
    return concat(covFuncOut.getParameters(), covFuncIn.getParameters());
}

void   SeparableGP::setParametersVector(vec p)
{
}

double SeparableGP::objective() const
{
    // TODO objective
    return 0.0;
}

vec    SeparableGP::gradient() const
{
    // TODO objective
    return zeros(1);
}

