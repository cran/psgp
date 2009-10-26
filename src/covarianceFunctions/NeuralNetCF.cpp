#include "NeuralNetCF.h"

/*
 * Constructor - pass in the length scale and variance
 */
NeuralNetCF::NeuralNetCF(double ls, double var, double offst)
    : CovarianceFunction("Neural network covariance function")
{
    // Make sure parameters are positive 
    assert( ls > 0.0 && var > 0.0 && offset > 0.0);  
    
    lengthScale = ls;   // Length scale
    variance = var;
    offset = offst;
    
    numberParameters = 3;
    setDefaultTransforms();
}


/*
 * Constructor - pass in a vector of parameters
 * (in this case [lengthscale, variance])
 */
NeuralNetCF::NeuralNetCF(vec parameters)
    : CovarianceFunction("Neural network covariance function")
{
    // Set number of parameters and check parameter vector has correct size
    numberParameters = 3;    
    assert(parameters.size() == getNumberParameters());

    // Set parameters
    assert(parameters(0) > 0.0 && parameters(1) > 0.0);
    lengthScale = parameters(0);
    variance = parameters(1);
    offset = parameters(2);

    setDefaultTransforms();
}


/**
 * Destructor
 */ 
NeuralNetCF::~NeuralNetCF()
{
}


/**
 * Return the name of the parameter of specified index
 */
string NeuralNetCF::getParameterName(int parameterNumber) const
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    switch (parameterNumber)
    {
    case 0: 
        return "Length scale";

    case 1:
        return "Variance";
        
    case 2:
        return "Offset";
    }
}


/**
 * Set given parameter to the specified value
 */
void NeuralNetCF::setParameter(int parameterNumber, const double value)
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
            
        case 2 :
            offset = value;
            break;
     }
}


/**
 * Return the parameter at index parameterNumber
 */
double NeuralNetCF::getParameter(int parameterNumber) const
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
        
        case 2 :
            return(offset);
            break;
    }

}


/**
 * Covariance between two points A and B
 */
inline double NeuralNetCF::computeElement(const vec& A, const vec& B) const
{
    if (A==B) 
        return computeDiagonalElement(A);
    
    vec Ahat = A/lengthScale;
    vec Bhat = B/lengthScale;
        
    double gamma = 1.0 / sqr(offset);
    
    return variance * asin( (gamma+dot(Ahat,Bhat)) / 
                            sqrt( (1.0+gamma+dot(Ahat,Ahat)) * (1.0+gamma+dot(Bhat,Bhat)) )
                          );
}

/**
 * Auto-covariance
 */
inline double NeuralNetCF::computeDiagonalElement(const vec& A) const
{
    double a = sum(sqr(A/lengthScale));
            
    double gamma = 1.0 / sqr(offset);
    
    return variance * asin( (gamma+a) / (1.0+gamma+a) );
}


/** 
 * Gradient of cov(X) w.r.t. given parameter number
 */
void NeuralNetCF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{
    assert(parameterNumber < getNumberParameters());
    assert(parameterNumber >= 0);

    Transform* t = getTransform(parameterNumber);
    double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

    switch(parameterNumber)
    {
        case 0 :
        {
            double gamma = 1.0 / sqr(offset);
            mat Xhat = X / lengthScale;
            
            vec Xi, Xj;
            for (int i=0; i<Xhat.rows(); i++)
            {
                Xi = Xhat.get_row(i);
                double vi = 1.0 + gamma + dot(Xi, Xi);
                
                for (int j=0; j<i; j++) 
                {
                    Xj = Xhat.get_row(j);
                    
                    double u  = gamma + dot(Xi,Xj);
                    double vj = 1.0 + gamma + dot(Xj, Xj);
                    // double vij = vi*vj;
                    // PD(i,j) = u * (vi+vj) / (vij * sqrt(vij - u*u) );
                    
                    double vij = sqrt(vi*vj);
                    PD(i,j) = (2*gamma*vij - (gamma+1)*u*(vi+vj)/vij) / (vij*sqrt(vij*vij-u*u));
                    PD(j,i) = PD(i,j);
                }
                // PD(i,i) = 2.0*(vi-1.0) / (vi*sqrt(2*vi-1.0));
                PD(i,i) = 2.0*(gamma+1.0-vi) / (vi*sqrt(2*vi-1.0));
            }
            PD *= variance * gradientModifier / lengthScale;
            break;
        }

        case 1 :
        {
            computeSymmetric(PD, X);
            PD *= (gradientModifier / variance);
            break;
        }

        case 2:
            double gamma = 1.0 / sqr(offset);
            mat Xhat = X / lengthScale;

            vec Xi, Xj;
            for (int i=0; i<Xhat.rows(); i++) 
            {
                Xi = Xhat.get_row(i);
                double vi = 1.0 + gamma + dot(Xi, Xi);

                for (int j=0; j<i; j++) 
                {
                    Xj = Xhat.get_row(j);

                    double u  = gamma + dot(Xi,Xj);
                    double vj = 1.0 + gamma + dot(Xj, Xj);
                    double vij = vi*vj;

                    PD(i,j) = (1.0 - 0.5*u*(vi+vj)/vij) / sqrt(vij-u*u); 
                    PD(j,i) = PD(i,j);
                }
                PD(i,i) = 1.0 / ( vi * sqrt(2*vi-1.0) );
            }
            PD *= - 2.0 * variance * gradientModifier / pow(offset,3.0);
            break;
    }

}
