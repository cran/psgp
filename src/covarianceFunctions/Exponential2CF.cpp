#include "Exponential2CF.h"

using namespace std;
using namespace itpp;

Exponential2CF::Exponential2CF(double l_x, double l_y, double l_xy, double var) : CovarianceFunction("Anisotropic Exponential")
{
	numberParameters = 4;
	setDefaultTransforms();
	setTransform(2, &id);

	L = zeros(2,2);
	L(0,0)  = l_x;
	L(1,1)  = l_xy;
	L(1,1)  = l_y;
	
	variance = var;
}


Exponential2CF::Exponential2CF(vec parameters) : CovarianceFunction("Anisotropic Exponential")
{
	numberParameters = 4;
	assert(parameters.size() == getNumberParameters());
	
	L = zeros(2,2);
	L(0,0)  = parameters(0);
	L(0,1)  = parameters(2);
	L(1,1)  = parameters(1);
	    
	variance = parameters(3);
	setDefaultTransforms();
	setTransform(2, &id);
}


Exponential2CF::~Exponential2CF()
{
}


inline double Exponential2CF::computeElement(const vec& A, const vec& B) const
{
    vec V = L*(A-B);
    return variance * exp( - 0.5 * sqrt( dot(V,V) ) );
}


inline double Exponential2CF::computeDiagonalElement(const vec& A) const
{
    return variance;
}


void Exponential2CF::setParameter(int parameterNumber, const double value)
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	
	switch(parameterNumber)
	{
	case 0: 
	    L(0,0) = value;
	    
	
	case 1: 
	    L(1,1) = value;
	    break;
	
	case 2:
	    L(0,1) = value;
	    break;
	    
	case 3: 
	    variance = value;
	    break;
	
	default: 
	    assert(false);
	    break;
	}
	
	
}

double Exponential2CF::getParameter(int parameterNumber) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
	case 0: 
	    return L(0,0);

	case 1:
	    return L(1,1);

	case 2:
	    return L(0,1);

	case 3 : 
	    return(variance);

	default: 
	    assert(false);
	    break;
	}
	cerr << "Warning: should not have reached here in GaussianCF::getParameter" << endl;
	return(0.0);
}


string Exponential2CF::getParameterName(int parameterNumber) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
	case 0 : 
	    return "Chol(Range)[0,0]";
	    
	case 1:
	    return "Chol(Range)[1,1]";
	    
	case 2:
	    return "Chol(Range)[0,1]";
	    
	case 3:
	    return "Variance";
	}
	return "Unknown parameter";
}

void Exponential2CF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	Transform* t = getTransform(parameterNumber);
	double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

	mat dL;
	
	switch(parameterNumber)
	{
	case 3:
	    computeSymmetric(PD, X);
        PD *= (gradientModifier / variance);
        return;
        
	case 0:
	    dL = zeros(2,2); dL(0,0) = 1.0;
	    break;
	
	case 1:
	    dL = zeros(2,2); dL(1,1) = 1.0;
	    break;
	
	case 2:
	    dL = zeros(2,2); dL(0,1) = 1.0;	
	    break;
	}
	
	
    vec v(2), dLv(2), Lv(2);
    double sLv2;
    
	for(int i=0; i<X.rows(); i++)
	{
	    for(int j=0; j<i; j++) 
	    {
	        v = X.get_row(i) - X.get_row(j);
	        dLv = dL*v;
	        Lv  = L*v;
	        sLv2 = sqrt( dot(Lv,Lv) );
	        PD(i,j) = -0.5 * dot(dLv, Lv) / sLv2 * exp( -0.5 * sLv2 ); 
	        PD(j,i) = PD(i,j); 
	    }
	    PD(i,i) = 0.0;
	}
	PD *= gradientModifier * variance;
}
