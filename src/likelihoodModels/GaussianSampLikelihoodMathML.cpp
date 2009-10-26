#include "GaussianSampLikelihoodMathML.h"

#include <cmath>


using namespace std;
using namespace itpp;
XERCES_CPP_NAMESPACE_USE

GaussianSampLikelihoodMathML::GaussianSampLikelihoodMathML(double Mean, double Variance, string MathMLString) : SamplingLikelihoodMathML(MathMLString)
{
	likelihoodMean = Mean;
	likelihoodVariance = Variance;
}

GaussianSampLikelihoodMathML::GaussianSampLikelihoodMathML(double Mean, double Variance, DOMDocument *doc) : SamplingLikelihoodMathML(doc)
{
	likelihoodMean = Mean;
	likelihoodVariance = Variance;
}

GaussianSampLikelihoodMathML::GaussianSampLikelihoodMathML(double Mean, double Variance, ExpressionNode *tree) : SamplingLikelihoodMathML(tree)
{
	likelihoodMean = Mean;
	likelihoodVariance = Variance;
}

GaussianSampLikelihoodMathML::~GaussianSampLikelihoodMathML()
{
}


double GaussianSampLikelihoodMathML::updateCoefficients(double& K1, double& K2, double Observation, double ModelMean, double ModelVariance) const
{

	double sigX = sqrt(ModelVariance);
	double y = Observation - likelihoodMean;
	
	// importance sampler mean and deviation
	double mq = ModelMean; // or (ModelMean + Observation) / 2
	double sq = 4 * sigX;   
	double ev = 0, pM = 0, pV = 0;
	
	// population monte carlo bit here
	for(int iCyc=0; iCyc < SamplingLikelihood::numberCycles; iCyc++)
	{
		vec sx = mq + (itpp::randn(SamplingLikelihood::numberSamples) * sq);
		// vec transH = SamplingLikelihood::modelFunction(sx);
		vec transH = apply_function(transform, sx);
		
		vec importanceWeights = (itpp::pow((sx - ModelMean) / sigX, 2.0) - itpp::pow((sx - mq) / sq, 2.0)) / 2.0;

		importanceWeights = itpp::exp(-importanceWeights) * sq / sigX;
		importanceWeights = itpp::elem_mult(importanceWeights, exp(-itpp::pow(y - transH, 2.0) / (2.0 * likelihoodVariance)) / sqrt(2.0 * itpp::pi * likelihoodVariance));
		ev = itpp::sum(importanceWeights);
		importanceWeights = importanceWeights / ev;
		pM = itpp::sum(itpp::elem_mult(sx, importanceWeights));		
		pV = itpp::sum(itpp::elem_mult(itpp::pow(sx, 2.0), importanceWeights)) - pow(pM, 2.0);
		mq = pM;
		sq = sqrt(abs(pV));
	}
	
	K1 = (pM - ModelMean) / ModelVariance;
	K2 = -(ModelVariance - pV) / pow(ModelVariance, 2.0);
	
	return log(ev);
	
}
