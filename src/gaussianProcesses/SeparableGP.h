#ifndef GP_H_
#define GP_H_

#include "ForwardModel.h"
#include "optimisation/Optimisable.h"
#include "covarianceFunctions/covarianceFunctions.h"

enum ParameterEstimationMethod { MARGINAL_LOG_LIKELIHOOD=0, MARGINAL_LOG_PRIOR=1 };

class SeparableGP : public ForwardModel, public Optimisable
{
public:
	SeparableGP(mat inputs, mat outputs, CovarianceFunction &cf_inputs, CovarianceFunction &cf_outputs);
	virtual ~SeparableGP();

	void makePredictions(vec &gpmean, mat &gpvar, mat new_inputs, CovarianceFunction& cf_inputs);
	void makePredictions(vec &gpmean, mat &gpvar, mat new_inputs);
	
	// Modifiers
	void setCovarianceFunctionInputs(CovarianceFunction &cf);
	void setCovarianceFunctionOutputs(CovarianceFunction &cf);
	void setParameterEstimationMethod(ParameterEstimationMethod method);
	
	// Accessors
	CovarianceFunction& getCovarianceFunctionInputs();
	CovarianceFunction& getCovarianceFunctionOutputs();
	string              getParameterEstimationMethodName();
	
	// Parameter estimation methods (Optimisable)
	vec    getParametersVector() const;
	void   setParametersVector(vec p);
	double objective() const;
	vec    gradient() const;
	
private:
    mat inputs;       // Inputs (training set)
    mat outputs;      // Outputs (training set)
    
    CovarianceFunction &covFuncIn;       // Covariance function between inputs
    CovarianceFunction &covFuncOut;      // Covariance function between outputs
    
    // mat *regressors(mat);               // Pointer to regressor (mean function)
    // mat beta;                           // Regression coefficients (mean function)
    
    

    // Parameter estimation functions
    ParameterEstimationMethod paramEstimationMethod;
    static string     parameterEstimationMethodName[];
    
    double marginalLogLikelihood(vec params);
    double marginalLogPosterior(vec params);
    
    vec    gradMarginalLogLikelihood(vec params);
    vec    gradMarginalLogPosterior(vec params);
    
    
};



// Initialise parameter estimation method names
string SeparableGP::parameterEstimationMethodName[] = { "Marginal log-likelihood", "Marginal log-prior" };



#endif /*GP_H_*/
