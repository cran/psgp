#include "Matern3CF.h"

/*
 * Constructor - pass in the length scale and variance
 */
Matern3CF::Matern3CF(double ls, double var)
    : CovarianceFunction("Matern 3/2 covariance function")
{
    // Make sure parameters are positive 
    assert( ls > 0.0 && var > 0.0 );  
    
    lengthScale = ls;
    variance = var; 
    
    numberParameters = 2;
    setDefaultTransforms();
}


/*
 * Constructor - pass in a vector of parameters
 * (in this case [lengthscale, variance])
 */
Matern3CF::Matern3CF(vec parameters)
    : CovarianceFunction("Matern 3/2 covariance function")
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
Matern3CF::~Matern3CF()
{
}


/**
 * Return the name of the parameter of specified index
 */
string Matern3CF::getParameterName(int parameterNumber) const
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
void Matern3CF::setParameter(int parameterNumber, const double value)
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
double Matern3CF::getParameter(int parameterNumber) const
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
inline double Matern3CF::computeElement(const vec& A, const vec& B) const
{
    if (A==B) 
        return computeDiagonalElement(A);
    
    double r = sqrt(3.0) * norm(A-B) / lengthScale;
     
    return variance * (1.0+r) * exp(-r);
}

/**
 * Auto-covariance
 */
inline double Matern3CF::computeDiagonalElement(const vec& A) const
{
    return variance;
}


/** 
 * Gradient of cov(X) w.r.t. given parameter number
 */
void Matern3CF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
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
            // computeSymmetric(PD, X);
            computeDistanceMatrix(R2, (sqrt(3.0) / lengthScale) * X);
            mat R = sqrt(R2); 
            elem_mult_out ( gradientModifier * (variance/lengthScale) * R2,  
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
