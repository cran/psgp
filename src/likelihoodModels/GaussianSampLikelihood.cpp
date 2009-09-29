#include "GaussianSampLikelihood.h"

#include <cmath>

using namespace std;
using namespace itpp;

GaussianSampLikelihood::GaussianSampLikelihood(double Mean, double Variance)
{
	likelihoodMean = Mean;
	likelihoodVariance = Variance;
}

GaussianSampLikelihood::~GaussianSampLikelihood()
{
}

// TODO  we need a variance check to make sure the variance doesn't collapse

double GaussianSampLikelihood::updateCoefficients(double& K1, double& K2, double Observation, double ModelMean, double ModelVariance) const
{

	double sigX = sqrt(ModelVariance);
	double y = Observation - likelihoodMean;

//	cout << "M: " << likelihoodMean << ", " << likelihoodVariance << ": " << numberSamples<< endl;
	
	// importance sampler mean and deviation
	double mq = ModelMean; // or (ModelMean + Observation) / 2
	

	// Inflate the prior in the initial population Monte Carlo run
	// This is an arbitrary factor.
	// sq =  4*sigX;
	double sq = 4 * sigX;   
	double ev = 0, pM = 0, pV = 0;
	
	// population monte carlo bit here
	
	//for nPop = 1:numPCycles;
	for(int iCyc=0; iCyc < SamplingLikelihood::numberCycles; iCyc++)
	{
		//sx = repmat(mq, nSample, 1) + randn(nSample,1 ) .* sq;
		//vec sx = mq + (samples * sq);
		vec sx = mq + (itpp::randn(SamplingLikelihood::numberSamples) * sq);
		
		//transH = feval(obsFunc, sx); 
//		vec transH = modelFunction(sx);
		vec transH = SamplingLikelihood::modelFunction(sx);
		
		//	computing weights
		//	impWeight = ((((sx-mX))./sigX).^2 - (((sx-mq))./sq).^2)/2;
		//  impWeight = exp(-impWeight)*sq/sigX;
		vec importanceWeights = (itpp::pow((sx - ModelMean) / sigX, 2.0) - itpp::pow((sx - mq) / sq, 2.0)) / 2.0;

		//cout << importanceWeights << endl;
		
		importanceWeights = itpp::exp(-importanceWeights) * sq / sigX;
		
		//cout << importanceWeights << endl;
		
		//	% multiply the importance weight with likelihood
		//	impWeight = impWeight.* exp(-(y - transH).^2/(2*covN))./(sqrt(2*pi*covN));
		//	ev        = sum(impWeight); % t1
		//	impWeight = impWeight./ev;
		importanceWeights = itpp::elem_mult(importanceWeights, exp(-itpp::pow(y - transH, 2.0) / (2.0 * likelihoodVariance)) / sqrt(2.0 * itpp::pi * likelihoodVariance));
		ev = itpp::sum(importanceWeights);

		
		importanceWeights = importanceWeights / ev;
		
		
		
		//	pM = sum(sx.*impWeight, 1); % t1/t3
		//	pV = sum((sx.^2).*impWeight) - pM.^2; % t2/t3 - (t1/t3)^2
		pM = itpp::sum(itpp::elem_mult(sx, importanceWeights));		
		pV = itpp::sum(itpp::elem_mult(itpp::pow(sx, 2.0), importanceWeights)) - pow(pM, 2.0);
		//	% Update the proposal distribution
		//	mq = pM;
		//	sq = sqrt(pV);
		mq = pM;
		sq = sqrt(abs(pV));
	}
	
//	%Sampling
//	K1S = (pM - mX) ./ sigX2;
//	K2S = -(sigX2-pV) ./ (sigX2).^2; 
//	logLikS = log(ev);

	K1 = (pM - ModelMean) / ModelVariance;
	K2 = -(ModelVariance - pV) / pow(ModelVariance, 2.0);
	
/*
	double sigX2 = ModelVariance + likelihoodVariance;
	double K2A = - 1 / sigX2;
	double K1A = -K2A * (Observation - modelFunction(ModelMean));
	double logLik = -(log(2 * M_PI * sigX2) + (Observation - modelFunction(ModelMean)) * K1A) / 2;
	//return logLik;
	*/
//	cout << "R: " << K1 << ", " << K2 << ", " << log(ev) << "  |  "  << K1A << ", " << K2A << ", " << logLik << endl;  

//	K1 = K1A;
//	K2 = K2A;
	
	return log(ev);
	
}
