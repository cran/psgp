#ifndef MATERN3CF_H_
#define MATERN3CF_H_

#include "CovarianceFunction.h"

#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>

/**
 * Isotropic Matern covariance function with nu=3/2.
 * 
 * (C) 2009 Remi Barillec <r.barillec@aston.ac.uk>
 */ 
class Matern3CF : public CovarianceFunction
{
public:
	Matern3CF(double lengthScale, double var);
	Matern3CF(vec parameters);
	virtual ~Matern3CF();
	
	// Compute covariance element cov(A,B) 
	inline double computeElement(const vec& A, const vec& B) const;
	
	// Compute autocovariance element cov(A,A)
    inline double computeDiagonalElement(const vec& A) const;
    
    // Partial derivative of covariance function with respect to a given parameter
    void    getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const;
    
    // Set and get a given parameter
    void    setParameter(int parameterNumber, const double value);
    double  getParameter(int parameterNumber) const;

	// Returns the name of a given parameter
    string  getParameterName(int parameterNumber) const;

private:
    static const double nu = 1.5;     // Nu parameter
    double lengthScale;                // Length scale parameter
    double variance;                   // Variance parameter
};


#endif /*MATERN3CF_H_*/
