#ifndef GENERICCF_H_
#define GENERICCF_H_

#include "CovarianceFunction.h"
#include "itppext/itppext.h"
#include <sstream>

using namespace itppext;

/**
 * Generic covariance function
 * 
 * Models a constant matrix of covariances, i.e. cov(X,X') = Sigma.
 * This is used to model the output covariance in the case of a 
 * separable multivariate Gaussian Process. 
 * 
 * The matrix is represented internally as a vector of lower triangular
 * elements from the Cholesky decomposition.  
 **/
class GenericCF : public CovarianceFunction
{
public:
	GenericCF(mat Sigma);
	virtual ~GenericCF();
	
	void   computeSymmetric(mat& C, const mat& X) const;
	void   computeSymmetricGrad(vec& V, const mat& X) const;
	void   computeCovariance(mat& C, const mat& X1, const mat& X2) const;
	void   computeDiagonal(mat& C, const mat& X) const;
	void   computeDiagonal(vec& C, const mat& X) const;   
	double computeElement(const vec& A, const vec& B) const;
	double computeDiagonalElement(const vec& A) const;

	void   getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const;

	void   setParameter(const int parameterNumber, const double value);
	double getParameter(const int parameterNumber) const;
	string getParameterName(const int parameterNumber) const;
	void   setParameters(const vec p);
	vec    getParameters() const;
	
private:
    int N;    // Number of columns/rows of Sigma
    vec Lvec; // Components of the Cholesky decomposition of Sigma
};

#endif /*GENERICCF_H_*/
