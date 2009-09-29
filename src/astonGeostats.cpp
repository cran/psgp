
#include <cassert>
#include <string>
#include <algorithm>

#include "astonGeostats.h"

#include "R.h"
#include "Rmath.h"


#include "itpp/itbase.h"
#include "itpp/itstat.h"

#include "gaussianProcesses/ForwardModel.h"
#include "gaussianProcesses/SequentialGP.h"
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
#include "covarianceFunctions/SumCovarianceFunction.h"
#include "covarianceFunctions/LogTransform.h"
#include "covarianceFunctions/IdentityTransform.h"
#include "covarianceFunctions/NegLogSigmoidTransform.h"
#include "covarianceFunctions/ConstantCF.h"

using namespace std;
using namespace itpp;



void learnParameters(int numObs, int numError, double *xData, double *yData, double *errorData, 
        int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, double *range, double *sill, 
        double *nugget, double *bias, int model)
{

	vec rndNums = itpp::randu(numObs);
	ivec rndIdx = itpp::sort_index(rndNums);

	int minLim = min(numObs, 400);

	mat X = itpp::zeros(minLim, 2);
	vec y = itpp::zeros(minLim);
	vec E = itpp::zeros(minLim);
	Vec<LikelihoodType *> noiseModels;

	// need to make this use a factory method from IT++
	// this would make things much more efficient
	for(int i=0; i < minLim; i++) {
		X(i, 0) = xData[rndIdx[i]];
		X(i, 1) = xData[rndIdx[i] + numObs];
		y(i) = yData[rndIdx[i]];
		E(i) = errorData[rndIdx[i]];		
	}

	// some reasonable starting parameters
	double r1=abs(min(max(X, 1) - min(X, 1)));
	double r2=abs(min(max(X, 2) - min(X, 2)));

	*range = 0.25 * ((r1+r2) / 2.0);
	*sill = abs(variance(y));
	*nugget = 0.5 * *sill;
	*bias = abs(1.0 / mean(y));

	// hardcode a gaussian covariance function so far
	ExponentialCF covComp1(*range, *sill);
  	WhiteNoiseCF covComp2(*nugget);
	ConstantCF covComp3(*bias);

	SumCovarianceFunction covSum(covComp1);
	covSum.addCovarianceFunction(covComp2);
	covSum.addCovarianceFunction(covComp3);

	cout << "===Starting Parameters===========================" <<endl;
	covSum.displayCovarianceParameters();

	// there will be a little memory leak here - should sort this out at somepoint...
	covSum.setTransform(0, new LogTransform());
	covSum.setTransform(1, new LogTransform());
	covSum.setTransform(2, new LogTransform());
	covSum.setTransform(3, new LogTransform());

	GaussianProcess gp(2, 1, X, y, covSum);

	SCGModelTrainer gpTrainer(gp);

	gpTrainer.setAnalyticGradients(true);
	gpTrainer.setCheckGradient(true);
	gpTrainer.Train(50);

	cout << "===Finishing Parameters===========================" <<endl;
	covSum.displayCovarianceParameters();

	*range  = covSum.getParameter(0);
	*sill   = covSum.getParameter(1);
	*nugget = covSum.getParameter(2);
	*bias   = covSum.getParameter(3);
	

}


void makePredictions(int numObs, int numPred, int numError, double *xData, double *yData,  double *xPred, 
							double *errorData, int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, 
							double *meanPred, double *varPred, double *range, double *sill, double *nugget,
							double *bias,int model)
{
    
	LikelihoodType *defaultLikelihood; 

	// leave until sort out the table
	defaultLikelihood = new GaussianLikelihood(*nugget * 0.01);

	// check the metadata table to see if there is anything available
	if(numMetadata == 0)
	{
		cout << "No noise models specified" << endl;
		cout << "Defaulting to GAUSSIAN:" << (*nugget * 0.01) << endl;
		defaultLikelihood = new GaussianLikelihood(*nugget * 0.01);
	}
	else
	{

		// parsing the metadata should be tidied up
		for(int i = 0; i < numMetadata; i++)
		{
			bool found = false;
			for(int j=0; j < numObs; j++)
			{
				if(errPtr[j] == (i + 1))
				{
					found = true;
					break;
				}
			} 

			if(found)
			{
				string s(metaDataTable[i]);

				string::size_type a = 0;  
				string::size_type b = s.find(',');
				string sensorName = s.substr(a, b-a);
				vec noiseParams;
				int idx = 0;
				a = ++b;  
				b = s.find(',', b);  

				while (b != string::npos)
				{
					noiseParams.set_size(idx + 1);
					noiseParams(idx++) = strtod((s.substr(a, b-a)).c_str(), NULL);
					a = ++b;  
					b = s.find(',', b);  
					if (b == string::npos)
					{
						noiseParams.set_size(idx + 1);
						noiseParams(idx++) = strtod((s.substr(a, s.length( ))).c_str(), NULL);
					}
				}  
			
				if(sensorName == "GAUSSIAN")
				{
					cout << (i + 1) << "gau:" << noiseParams(1) << "; " << endl;
				}
				else
				{
					cout << "Unrecognised sensor model (" << sensorName << ")- defaulting to Gaussian" << endl;
				}
			}		
		}
	}

	vec rndNums = itpp::randu(numObs);
	ivec rndIdx = itpp::sort_index(rndNums);

	int minLim = min(numObs, 500);

	mat Xpred = itpp::zeros(numPred, 2);
	vec predY = itpp::zeros(numPred);
	vec predVar = itpp::zeros(numPred);

	mat X = itpp::zeros(minLim,2);
	vec y = itpp::zeros(minLim);
	vec E = itpp::zeros(minLim);
	
	// need to make this use a factory method from IT++
	// this would make things much more efficient
	for(int i=0; i < minLim; i++) {
		X(i, 0) = xData[rndIdx[i]];
		X(i, 1) = xData[rndIdx[i] + numObs];
		y(i) = yData[rndIdx[i]];
		E(i) = errorData[rndIdx[i]];
	}

	
	for(int i=0; i < numPred; i++) {
		Xpred(i, 0) = xPred[i];
		Xpred(i, 1) = xPred[i + numPred];
	}


	// hardcode a gaussian covariance function so far
	ExponentialCF covComp1(*range, *sill);
  	WhiteNoiseCF covComp2(*nugget);
	ConstantCF covComp3(*bias);

	SumCovarianceFunction covSum(covComp1);
	covSum.addCovarianceFunction(covComp2);
	covSum.addCovarianceFunction(covComp3);

	SumCovarianceFunction covSumPred(covComp1);
	covSumPred.addCovarianceFunction(covComp3);
	
	// Set number of active points to max(400, size(obs))
	int n_active = max(100, numObs);
	SequentialGP ssgp(2,1,n_active,X,y,covSum,1);

	cout << "Computing posterior..." << endl;
	ssgp.computePosterior(*defaultLikelihood);

	cout << "Predicting..."<<endl;

	// should add an option here to select one or the other...
	if(0)
	{
	    ssgp.makePredictions(predY, predVar, Xpred, covSumPred);
	}
	else
	{
	    int startVal = 0;
	    int chunkSize = 1000;
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

	        // predYChunk.set_size(chunkSize);
	        // predVarChunk.set_size(chunkSize);
	        vec predYChunk(chunkSize);
	        vec predVarChunk(chunkSize);

	        ssgp.makePredictions(predYChunk, predVarChunk, XpredChunk, covSumPred);

	        predY.replace_mid(startVal, predYChunk);
	        predVar.replace_mid(startVal, predVarChunk);
	        
	        /* RB: Same as above - used for debugging 
	        for (int i=0; i<chunkSize; i++) 
	        {
	            predY(startVal+i) = predYChunk(i);
	            predVar(startVal+i) = predVarChunk(i);
	        }
	        */
	        
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

	covSum.displayCovarianceParameters();

}

