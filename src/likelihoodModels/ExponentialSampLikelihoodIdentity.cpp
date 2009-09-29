#include "ExponentialSampLikelihoodIdentity.h"

#include <cmath>

using namespace std;
using namespace itpp;

ExponentialSampLikelihoodIdentity::ExponentialSampLikelihoodIdentity(double Lambda) : ExponentialSampLikelihood(Lambda)
{
}

ExponentialSampLikelihoodIdentity::~ExponentialSampLikelihoodIdentity()
{
}

double ExponentialSampLikelihoodIdentity::modelFunction(const double x) const
{
// identity bias
	return(x);
}
