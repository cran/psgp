#include "Matern5CF.h"

/*
 * Constructor - pass in the length scale
 */
Matern5CF::Matern5CF(double ls, double var)
    : CovarianceFunction("Matern 5/2 covariance function")
{
    // Make sure length scale is positive 
    assert( ls > 0.0 && var > 0.0 );  
    
    lengthScale = ls;   // Length scale
    variance = var;
    
    numberParameters = 2;
    setDefaultTransforms();
}


/*
 * Constructor - pass in a vector of parameters
 * (in this case [lengthscale, variance])
 */
Matern5CF::Matern5CF(vec parameters)
    : CovarianceFunction("Matern 5/2 covariance function")
{
    // Set number of parameters and check parameter vector has correct size
    numberParameters = 2;    
    assert(parameters.size() == getNumberParameters());

    // Set parameters
    assert(parameters(0) > 0.0 && parameters(1) > 0.0);
    lengthScale = parameters(0);
    variance = parameters(1);

    setDefaultTransforms();
}


/*
 * Destructor
 */
Matern5CF::~Matern5CF()
{
}


/**
 * Return the name of the parameter of specified index
 */
string Matern5CF::getParameterName(int parameterNumber) const
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    switch (parameterNumber)
    {
    case 0: 
        return "Length scale";

    case 1:
        return "Variance";
    }
}


/**
 * Set given parameter to the specified value
 */
void Matern5CF::setParameter(int parameterNumber, const double value)
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    switch(parameterNumber)
    {
        case 0 : 
            lengthScale = value;
            break;
            
        case 1 : 
            variance = value;
            break;
     }
}


/**
 * Return the parameter at index parameterNumber
 */
double Matern5CF::getParameter(int parameterNumber) const
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    switch(parameterNumber)
    {
        case 0 : 
            return(lengthScale);
            break;
            
        case 1 : 
            return(variance);
            break;
    }

}


/**
 * Covariance between two points A and B
 */
inline double Matern5CF::computeElement(const vec& A, const vec& B) const
{
    if (A==B) 
        return computeDiagonalElement(A);
    
    double r = sqrt(5.0) * norm(A-B) / lengthScale;
     
    return variance * ( 1.0 + r + sqr(r)/3.0 ) * exp(-r);
}

/**
 * Auto-covariance
 */
inline double Matern5CF::computeDiagonalElement(const vec& A) const
{
    return variance;
}


/** 
 * Gradient of cov(X) w.r.t. given parameter number
 */
void Matern5CF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    Transform* t = getTransform(parameterNumber);
    double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

    switch(parameterNumber)
    {
        case 0 :
        {
            mat R2(PD.rows(), PD.cols());
            computeDistanceMatrix(R2, (sqrt(5.0) / lengthScale) * X);
            mat R = sqrt(R2); 
            elem_mult_out ( gradientModifier * (variance/(3.0*lengthScale)) * elem_mult(R2, 1.0+R),
                            exp(-R), 
                            PD
            );  
            
            break;
        }

        case 1 :
        {
            computeSymmetric(PD, X);
            PD *= (gradientModifier / variance);
            break;
        }
    }

}

