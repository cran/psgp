#include "GaussianLikelihood.h"

#include <cmath>

using namespace std;

GaussianLikelihood::GaussianLikelihood(const double lp)
{
	 likelihoodParameter = lp;
}

GaussianLikelihood::~GaussianLikelihood()
{
}


double GaussianLikelihood::updateCoefficients(double& K1, double& K2, const double Observation, 
											  const double ModelMean, const double ModelVariance) const
{
	// Lehel: Sec. 2.4
	double sigX2 = ModelVariance + likelihoodParameter;
	double logLik;
	K2 = - 1 / sigX2;
	K1 = -K2 * (Observation - ModelMean);
	logLik = -(log(2 * M_PI * sigX2) + (Observation - ModelMean) * K1) / 2;
	return logLik;
}
