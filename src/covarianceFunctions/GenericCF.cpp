#include "GenericCF.h"

/**
 * Constructor
 * @param Sigma The covariance matrix)
 */
GenericCF::GenericCF(mat Sigma)
: CovarianceFunction("Generic covariance function (fixed covariance matrix)")
{
    // Store the coefficients of the Cholesky decomposition of Sigma
    Lvec = ltr_vec(chol(Sigma)); 

    numberParameters = Lvec.length();
    setDefaultTransforms();
}

/**
 * Destructor
 */
GenericCF::~GenericCF()
{
}

/**
 * Sets a given parameter to the specified value
 */
void   GenericCF::setParameter(const int parameterNumber, const double value)
{
    Lvec(parameterNumber) = value;
}

/**
 * Gets parameter by index
 * 
 * @param parameterNumber Index of the parameter
 * @return The parameter with index parameterNumber
 */
double GenericCF::getParameter(const int parameterNumber) const
{
    return Lvec(parameterNumber);    
}

/**
 * Gets the name of parameter at specified index
 * 
 * @returns The name of the parameter
 **/ 
string GenericCF::getParameterName(const int parameterNumber) const 
{
    stringstream s;
    s << "Parameter " << parameterNumber << " in Cholesky decomposition of covariance matrix";
    
    return s.str();
}

/**
 * Sets all parameters from vector
 */
void GenericCF::setParameters(const vec p)
{
    Lvec = p;
}

/**
 * Returns a vector of parameters
 */
vec GenericCF::getParameters() const
{
    return Lvec;
}





/**
 * Computes the auto-covariance matrix of inputs X. 
 * In this case, this simply returns Sigma (X is here for compatibility).
 */
void   GenericCF::computeSymmetric(mat& Sigma, const mat& X) const 
{
    mat L = utr_mat(Lvec);          // Cholesky decomposition of Sigma
    Sigma = L*L.transpose();
}

void   GenericCF::computeCovariance(mat& Sigma, const mat& X1, const mat& X2) const 
{
    cout << "**Warning: GenericCF::computeCovariance" << endl;
    cout << "**         You should not be using this method." << endl;
    
    computeSymmetric(Sigma, X1);
}


void   GenericCF::computeSymmetricGrad(vec& V, const mat& X) const {}

void   GenericCF::computeDiagonal(mat& Sigma, const mat& X) const 
{
    vec d;
    computeDiagonal(d, X);     // Compute vector of diagonal elements
    Sigma = diag(d);           // And assign it to diagonal of Sigma
}

void   GenericCF::computeDiagonal(vec& Sigma, const mat& X) const 
{
    Sigma = zeros(N);
    int k = 0;
    for (int i=0; i<N; i++) {
        Sigma(i) = dot(Lvec(k,k+i), Lvec(k,k+i));
    }
          
}

double GenericCF::computeElement(const vec& A, const vec& B) const 
{ 
    // TODO computeElement
    return 0.0; 
}

double GenericCF::computeDiagonalElement(const vec& A) const 
{ 
    // TODO computeDiagonalElement
    return 0.0; 
}

/**
 * Returns the gradient of Sigma with respect to parameter at specified index
 */
void   GenericCF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const 
{
    vec dLvec = zeros(Lvec.length());
    dLvec(parameterNumber) = 1.0;
    
    mat L = ltr_mat(Lvec);
    mat dL = ltr_mat(dLvec);
    PD = dL*L.transpose() + L*dL.transpose();
}