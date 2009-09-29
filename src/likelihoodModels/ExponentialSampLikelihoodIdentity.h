#ifndef EXPONENTIALSAMPLIKELIHOODIDENTITY_H_
#define EXPONENTIALSAMPLIKELIHOODIDENTITY_H_

#include "LikelihoodType.h"
#include "ExponentialSampLikelihood.h"


class ExponentialSampLikelihoodIdentity : public ExponentialSampLikelihood
{
public:
	ExponentialSampLikelihoodIdentity(double Lambda);
	virtual ~ExponentialSampLikelihoodIdentity();
	
	double modelFunction(const double x) const;

};

#endif /*EXPONENTIALSAMPLIKELIHOODIDENTITY_H_*/
