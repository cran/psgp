#include "MaternCF.h"

using namespace std;
using namespace itpp;

MaternCF::MaternCF(double lengthscale, double var, double kappa) : CovarianceFunction("Isotropic Matern")
{
	numberParameters = 3;
	setDefaultTransforms();

	range = lengthscale;
	variance = var;
	smoothness = kappa;
}

MaternCF::MaternCF(vec parameters) : CovarianceFunction("Isotropic Matern")
{
	numberParameters = 3;
	assert(parameters.size() == getNumberParameters());
	range = parameters(0);
	variance = parameters(1);
	smoothness = parameters(2);
	setDefaultTransforms();
}


MaternCF::~MaternCF()
{
}

inline double MaternCF::computeElement(const vec& A, const vec& B) const
{
	return calcMatern(A - B);
}

inline double MaternCF::computeDiagonalElement(const vec& A) const
{
	return calcMaternDiag();
}

inline double MaternCF::calcMatern(const vec& V) const
{
	double covf;
	double distMat = sqrt(sum(pow(V, 2.0))) / (2*range);
	double dist = sqrt(2 * smoothness) * distMat;
	
	double logCov = log(computeBessel(dist, smoothness));

	logCov -= dist;
	logCov += smoothness * log(dist);
	logCov += log(variance);
	logCov -= (smoothness - 1) * log(2);
	logCov -= log(itpp::gamma(smoothness));
	covf = exp(logCov);

	if(isnan(covf))
	{
		cout << "Negative covariance" << endl;
		covf = variance;
	}
	else
	{
		if(isinf(covf)) {
			cout << "Infinite covariance" << endl;
			covf = variance;
		}
	}

	return covf;
}

inline double MaternCF::calcGaussian(const vec& V) const
{
	return variance * exp( -0.5 * sum(itpp::pow(V , 2.0) / pow(range, 2.0)));
}

inline double MaternCF::calcMaternDiag() const
{
	return variance;
}

void MaternCF::setParameter(int parameterNumber, const double value)
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
		case 0 : range = value;
					break;
		case 1 : variance = value;
					break;
		case 2 :	if(smoothness > 100.0)
					{
						smoothness = 100.0; // constrain smoothness parameter
					}
					smoothness = value;
					break;
		default: assert(false);
					break;
	}
}

double MaternCF::getParameter(int parameterNumber) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
		case 0 : return(range);
					break;
		case 1 : return(variance);
					break;
		case 2 : return(smoothness);
					break;
		default: assert(false);
					break;
	}
	cerr << "Warning: should not have reached here in GaussianCF::getParameter" << endl;
	return(0.0);
}

string MaternCF::getParameterName(int parameterNumber) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	switch(parameterNumber)
	{
		case 0 : return("Range");
					break;
		case 1 : return("Variance");
					break;
		case 2 : return("Smoothness");
					break;
		default: break;

	}
	return("Unknown parameter");
}

/*
if nu<100;
  iS = find(ddd<1e-5);
  iL = setdiff([1:length(ddd)],iS);
  % computing for large distances
  dPlus(iL,:) = - 2*nu./ddd(iL) ...
      + besselk(nu+1,ddd(iL),1)./besselk(nu, ddd(iL), 1);
  dPlus(iS,:) = 0;
else
  dPlus = ddd./2./nu;
end;
*/
void MaternCF::getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const
{
	assert(parameterNumber < getNumberParameters());
	assert(parameterNumber >= 0);

	Transform* t = getTransform(parameterNumber);
	double gradientModifier = t->gradientTransform(getParameter(parameterNumber));

	switch(parameterNumber)
	{
		case 0 :
		{
			mat DM(PD.rows(), PD.cols());
			computeSymmetric(PD, X);
//			computeDistanceMatrix(DM, X);
			DM  = PD / variance;
//			DM = sqrt(DM  * ((2 * smoothness) / pow(range, 2.0)));
//			cout << DM << endl;
//			elem_mult_inplace(DM, PD);
			mat res(PD.rows(), PD.cols());
			
			if(smoothness < 100)
			{
			
				for(int i = 0; i < X.rows(); i++)
				{
					for(int j = 0; j < X.cols(); j++)
					{
						if(DM(i, j) < (1e-5))
						{
							res(i, j) = 0.0;
						}
						else
						{
							double val1 = computeBessel(DM(i, j), smoothness + 1);
							double val2 = computeBessel(DM(i, j), smoothness);
							res(i, j) = ((-2 * smoothness) / DM(i, j)) + (val1 / val2);
						}
					}
				}			
			}
			else
			{
				res = DM / (2 * smoothness);
			}
			
			elem_mult_inplace(res, PD);
			
			return;
			break;
		}

		case 1 :
		{
			computeSymmetric(PD, X);
			PD *= (gradientModifier / variance);
			return;
			break;
		}
		case 2 :
		{
			// we cannot compute gradient for smoothness parameter
			PD = zeros(PD.rows(), PD.cols());
			return;	
			break;
		}
		default: break;
	}
	cerr << "Warning: should not have reached here in MaternCF::getParameterPartialDerivative" << endl;
}

inline double MaternCF::computeBessel(const double x, const double smoothness) const
{
	#ifdef USING_R
	return bessel_k(x, smoothness, 1.0);
	#else
	return gsl_sf_bessel_Knu_scaled(smoothness, x);
	#endif
}



mat MaternCF::besselkMat(const mat& X, double nu) const
{
	mat PD(X.rows(), X.cols());
	for(int i = 0; i < X.rows(); i++)
	{
		for(int j = 0; j < X.cols(); j++)
		{
			PD(i, j) = computeBessel(X(i, j), nu);
		}
	}
	return PD;
}




