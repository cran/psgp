#include "CovarianceFunction.h"

using namespace std;
using namespace itpp;

CovarianceFunction::CovarianceFunction(string name)
{
	covarianceName = name;
	transformsApplied = false;
}

CovarianceFunction::~CovarianceFunction()
{

}


/**
 * Compute auto-covariance of a single input X  
 */ 
void CovarianceFunction::computeSymmetric(double &c, const vec& x) const 
{
   mat C(1,1);
   computeSymmetric(C,x.transpose());
   c = C(0,0);
}


void CovarianceFunction::computeSymmetric(mat& C, const mat& X) const
{
	// ensure that data dimensions match supplied covariance matrix
	assert(C.rows() == X.rows());
	assert(C.cols() == X.rows());

	if (X.rows() == 1)
	{
	    C.set(0, 0, computeDiagonalElement(X.get_row(0)));
	    return;
	}
	
	// calculate the lower and upper triangles
	double d;
	
	for(int i=0; i<X.rows() ; i++)
	{
		for(int j=0; j<i; j++)
		{
		    d = computeElement(X.get_row(i), X.get_row(j));
		    C.set(i, j, d);
			C.set(j, i, d);
		}
	}

	// calculate the diagonal part
	for(int i=0; i<X.rows() ; i++)
	{
		C.set(i, i, computeDiagonalElement(X.get_row(i)));
	}
}

void CovarianceFunction::computeSymmetricGrad(vec& V, const mat& X) const
{

}



/**
 * Compute covariance between a set of inputs X and single input x
 */
void CovarianceFunction::computeCovariance(vec& c, const mat& X, const vec& x) const
{
    mat Xmat(x);
    mat C(X.rows(),1);
    computeCovariance(C,X,Xmat.transpose());
    c = C.get_col(0);
}

/**
 *  Covariance between points (rows) in X1 and points in X2
 */
void CovarianceFunction::computeCovariance(mat& C, const mat& X1, const mat& X2) const
{
	assert(C.rows() == X1.rows());
	assert(C.cols() == X2.rows());

	for(int i=0; i<X1.rows() ; i++)
	{
		for(int j=0; j<X2.rows(); j++)
		{
			C.set(i, j, computeElement(X1.get_row(i), X2.get_row(j)));
		}
	}
}

void CovarianceFunction::computeDiagonal(mat& C, const mat& X) const
{
	// calculate the diagonal part
	for(int i=0; i<X.rows() ; i++)
	{
		C.set(i, i, computeDiagonalElement(X.get_row(i)));
	}
}

void CovarianceFunction::computeDiagonal(vec& C, const mat& X) const
{
	// calculate the diagonal part
	for(int i=0; i<X.rows() ; i++)
	{
		C.set(i, computeDiagonalElement(X.get_row(i)));
	}
}

/**
 * Display information about current covariance function
 * If wanted, the output can be indented by a number space characters
 * as specified in argument.
 */
void CovarianceFunction::displayCovarianceParameters(int nspaces) const
{
	cout.setf(ios::fixed);
	cout.precision(4);

	string space = string(nspaces, ' ');
	
	cout << space << "Covariance function : " << covarianceName << endl;

	vec t = getParameters();

	for(int i=0; i < (t.size()); i++)
	{
		cout << space << getParameterName(i) << " (P" << (i) << ") : ";
		cout << (transforms[i]->backwardTransform(t(i)));
		if(transformsApplied)
		{
			cout  << " (" << transforms[i]->type() << " transformed)";
		}
		else
		{
			cout << " (no transform applied)";
		}
		cout << endl;
	}
	// cout << "====================" << endl;

}

void CovarianceFunction::setParameters(const vec p)
{
	assert(transforms.size() == p.size());
	for(int i = 0; i < getNumberParameters() ; i++)
	{
		setParameter(i, transforms[i]->backwardTransform(p(i)));
	}
}

vec CovarianceFunction::getParameters() const
{
	assert(transforms.size() == numberParameters);
	vec result;
	result.set_size(getNumberParameters());
	for(int i = 0; i < getNumberParameters() ; i++)
	{
		result[i] = transforms[i]->forwardTransform(getParameter(i));
	}
	return result;
}

int CovarianceFunction::getNumberParameters() const
{
	return numberParameters;
}

void CovarianceFunction::setTransform(int parameterNumber, Transform* newTransform)
{
	assert(parameterNumber >= 0);
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber < transforms.size());
	transforms[parameterNumber] = newTransform;
}

void CovarianceFunction::setDefaultTransforms()
{
	int sz = getNumberParameters();
	LogTransform* defaultTransform = new LogTransform();
	// IdentityTransform* defaultTransform = new IdentityTransform();

	for(int i = transforms.size(); i < sz; i++)
	{
		transforms.push_back(defaultTransform);
	}
	transformsApplied = true;
}

Transform* CovarianceFunction::getTransform(int parameterNumber) const
{
	assert(parameterNumber >= 0);
	assert(parameterNumber < getNumberParameters());
	return transforms[parameterNumber];
}

void CovarianceFunction::computeDistanceMatrix(mat& DM, const mat& X) const
{
	double value;
	assert(DM.rows() == X.rows());
	assert(DM.cols() == X.rows());

	for(int i=0; i<X.rows() ; i++)
	{
		for(int j=0; j<i; j++)
		{
			value = sum( sqr( X.get_row(i) - X.get_row(j) ) );
			DM.set(i, j, value);
			DM.set(j, i, value);
		}
		DM.set(i, i, 0.0);
	}
}
