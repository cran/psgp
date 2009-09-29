#include "ConstantCF.h"

ConstantCF::ConstantCF(double amp) : CovarianceFunction("Constant")
{
	numberParameters = 1;
	setDefaultTransforms();
	amplitude = amp;

}

ConstantCF::~ConstantCF()
{

}

inline double ConstantCF::computeElement(const vec& A, const vec& B) const
{
	return 1 / amplitude;
}

inline double ConstantCF::computeDiagonalElement(const vec& A) const
{
	return 1 / amplitude;
}

double ConstantCF::getParameter(int parameterNumber) const
{
	assert(parameterNumber == 0);

	switch(parameterNumber)
	{
		case 0 : return(amplitude);
					break;
		default: break;
	}
	cerr << "Warning: should not have reached here in ConstantCF::getParameter" << endl;
	return(0.0);
}

void ConstantCF::setParameter(int parameterNumber, const double value)
{
	assert(parameterNumber == 0);

	switch(parameterNumber)
	{
		case 0 : amplitude = value;
					break;
		default: break;
	}
}

string ConstantCF::getParameterName(int parameterNumber) const
{
	assert(parameterNumber == 0);

	switch(parameterNumber)
	{
		case 0 : return("Amplitude");
					break;
		default: break;

	}
	return("Unknown parameter");
}

void ConstantCF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{
	assert(parameterNumber == 0);

	Transform* t = getTransform(parameterNumber);
	double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

	switch(parameterNumber)
	{
		case 0 :
		{
			for(int i = 0; i < X.rows(); i++)
			{
				for(int j = 0; j < X.rows(); j++)
				{
					PD(i,j) = -1 / amplitude;
				}
			}
			return;
			break;
		}
	}
	cerr << "Warning: should not have reached here in ConstantCF::getParameterPartialDerivative" << endl;
}
