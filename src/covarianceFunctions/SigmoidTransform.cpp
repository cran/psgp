#include "SigmoidTransform.h"

SigmoidTransform::SigmoidTransform()
{
	transformName = "Sigmoid";
}

SigmoidTransform::~SigmoidTransform()
{
}

double SigmoidTransform::forwardTransform(const double a) const
{
	//cout << "NegLogSigmoidTransform::forwardTransform" << endl;
	return log(a / (1.0 - a));	// inv sigmoid

}

double SigmoidTransform::backwardTransform(const double b) const
{
	//cout << "NegLogSigmoidTransform::backwardTransform" << endl;
	return (1.0 / (1.0 + exp(-b)));	// sigmoid
}

double SigmoidTransform::gradientTransform(const double g) const
{
	return (g * (1 - g));
}
