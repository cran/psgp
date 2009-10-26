#include "SumCovarianceFunction.h"

using namespace std;
using namespace itpp;

SumCovarianceFunction::SumCovarianceFunction(vector<CovarianceFunction> cfVec) : CovarianceFunction("Sum Covariance")
{

	cout << "NOT IMPLEMENTED YET!!!" << endl;

}

SumCovarianceFunction::SumCovarianceFunction(CovarianceFunction& cf) : CovarianceFunction("Sum Covariance")
{
	covFunctions.push_back(&cf);
	numberParameters = cf.getNumberParameters();
	setDefaultTransforms();
}

void SumCovarianceFunction::addCovarianceFunction(CovarianceFunction& cf)
{
	covFunctions.push_back(&cf);
	numberParameters = numberParameters + cf.getNumberParameters();
	setDefaultTransforms();
}

SumCovarianceFunction::~SumCovarianceFunction()
{
}

inline double SumCovarianceFunction::computeElement(const vec& A, const vec& B) const
{
	double k = 0.0;

	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		k = k + covFunctions[i]->computeElement(A, B);
	}

	return k;
}

inline double SumCovarianceFunction::computeDiagonalElement(const vec& A) const
{
	double k = 0.0;

	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		k = k + covFunctions[i]->computeDiagonalElement(A);
	}

	return k;
}

void SumCovarianceFunction::displayCovarianceParameters(int nspaces) const
{
	string space(nspaces, ' ');
	
    cout << space << "Covariance function : Sum" << endl;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		cout << space << "+ Component: " << (i+1) << endl;
		covFunctions[i]->displayCovarianceParameters(nspaces+2);
	}
}

void SumCovarianceFunction::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{

//	vec result;
	int pos = 0;

//	result.set_size(getNumberParameters());

	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		for(int j = 0; j < (covFunctions[i]->getNumberParameters()) ; j++)
		{
			if(parameterNumber == pos)
			{
				covFunctions[i]->getParameterPartialDerivative(PD, j, X);
				return;
			}
			pos = pos + 1;
		}
	}
}


Transform* SumCovarianceFunction::getTransform(int parameterNumber) const
{
	int pos = 0;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		for(int j = 0; j < (covFunctions[i]->getNumberParameters()) ; j++)
		{
			if(parameterNumber == pos)
			{
				return covFunctions[i]->getTransform(j);
			}
			pos = pos + 1;
		}
	}

}


void SumCovarianceFunction::setTransform(int parameterNumber, Transform* newTransform)
{
	assert(parameterNumber >= 0);
	assert(parameterNumber < getNumberParameters());

	//transforms[parameterNumber] = newTransform;

	int pos = 0;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		for(int j = 0; j < (covFunctions[i]->getNumberParameters()) ; j++)
		{
			if(parameterNumber == pos)
			{
				covFunctions[i]->setTransform(j, newTransform);
				return;
			}
			pos = pos + 1;
		}
	}
}


void SumCovarianceFunction::setParameters(const vec p)
{
	int parFrom = 0;
	int parTo = 0;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
	    // RB: Extract parameters for covariance function j
	    parFrom = parTo;
	    parTo += covFunctions[i]->getNumberParameters();
	    covFunctions[i]->setParameters( p(parFrom, parTo-1) ); 
	        
	    /*
	    for(int j = 0; j < (covFunctions[i]->getNumberParameters()) ; j++)
		{
			Transform* t = covFunctions[i]->getTransform(j);
			double d = t->backwardTransform(p(pos));
			covFunctions[i]->setParameter(j, d);
			pos = pos + 1;
		}
	    */
	}
}

vec SumCovarianceFunction::getParameters() const
{
	vec result;
	int pos = 0;

	result.set_size(getNumberParameters());

	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		for(int j = 0; j < (covFunctions[i]->getNumberParameters()) ; j++)
		{
			Transform* t = covFunctions[i]->getTransform(j);
			double d = t->forwardTransform(covFunctions[i]->getParameter(j));
			result[pos] = d;
			pos = pos + 1;
		}
	}
	return result;
}

void SumCovarianceFunction::setParameter(const int parameterNumber, const double value)
{

	cout << "SumCovarianceFunction::setParameter" << endl;
/*
	int pnum = 0;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		int numParams = covFunctions[i]->getNumberParameters();
		if(parameterNumber < (pnum + numParams))
		{
			covFunctions[i]->setParameter(parameterNumber - pnum, value);
			return;
		}
		pnum = pnum + numParams;
	}
	*/
	cerr << "We shouldn't reach here - setParam : " << parameterNumber << endl;
}

double SumCovarianceFunction::getParameter(const int parameterNumber) const
{
	int pos = 0;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		for(int j = 0; j < (covFunctions[i]->getNumberParameters()) ; j++)
		{
			if(parameterNumber == pos)
			{
				return covFunctions[i]->getParameter(j);
			}
			pos = pos + 1;
		}
	}
	assert(false);
	return(0.0);

}

string SumCovarianceFunction::getParameterName(const int parameterNumber) const
{
	int pnum = 0;
	for(std::vector<CovarianceFunction *>::size_type i = 0; i < covFunctions.size(); i++)
	{
		int numParams = covFunctions[i]->getNumberParameters();
		if(parameterNumber < (pnum + numParams))
		{
			return (covFunctions[i]->getParameterName(parameterNumber - pnum));
		}
		pnum = pnum + numParams;
	}
	cerr << "We shouldn't reach here - getParamName" << endl;
	return("Unknown");
}

