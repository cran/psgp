#include "NegLogSigmoidTransform.h"

NegLogSigmoidTransform::NegLogSigmoidTransform()
{
	transformName = "Negative Log Sigmoid";
}

NegLogSigmoidTransform::~NegLogSigmoidTransform()
{
}

double NegLogSigmoidTransform::forwardTransform(const double a) const
{
	if(a < -MAXEXP)
	{
		return log(exp(a) - 1.0);
	}
	return log(exp(a) - 1.0);
}

double NegLogSigmoidTransform::backwardTransform(const double b) const
{

	if(b < -MAXEXP)
	{
		return eps;
	}
	else
	{
		if(b > MAXEXP)
		{
			return b;
		}
	}

	return log(1.0 + exp(b));
}

double NegLogSigmoidTransform::gradientTransform(const double g) const
{

	if(g > MAXEXP)
	{
		return 1.0;
	}
	return ((exp(g) - 1.0) / exp(g));
}
