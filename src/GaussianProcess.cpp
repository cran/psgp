#include "GaussianProcess.h"

using namespace std;
using namespace itpp;

GaussianProcess::GaussianProcess(int Inputs, int Outputs, mat& Xdata, vec& ydata, CovarianceFunction& cf) : Locations(Xdata), Observations(ydata), covFunc(cf), ForwardModel(Inputs, Outputs)
{
	assert(Locations.rows() == Observations.size());

}

GaussianProcess::~GaussianProcess()
{
}

//void GaussianProcess::train(mat& C, const mat& X) const
//{
//	mat A;
//	A.set_size(C.rows(), C.cols(), false);
//	covFunc->computeSymmetric(A, X);
//	inv(A, C);
//}

void GaussianProcess::makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const
{
	assert(Mean.size() == Variance.size());
	assert(Xpred.rows() == Mean.size());

	mat Sigma(Observations.size(), Observations.size());
	mat Cpred(Locations.rows(), Xpred.rows());

	covFunc.computeCovariance(Cpred, Locations, Xpred);
	covFunc.computeSymmetric(Sigma, Locations);

	vec alpha = ls_solve(Sigma, Observations);
	Mean = Cpred.transpose() * alpha;

	vec variancePred(Xpred.rows());
	covFunc.computeDiagonal(variancePred, Xpred);
	mat v = ls_solve(computeCholesky(Sigma), Cpred);
//	cout << sum(elem_mult(v, v)) << endl;
	Variance = variancePred - sum(elem_mult(v, v));
	

}

void GaussianProcess::makePredictions(vec& Mean, vec& Variance, const mat& Xpred, const mat& C) const
{

	assert(Mean.size() == Variance.size());
	assert(Xpred.rows() == Mean.size());

	mat Sigma, cholSigma;
	cholSigma.set_size(Observations.size(), Observations.size());
	
	Sigma.set_size(Observations.size(), Observations.size());
//	cholSigma.set_size(Observations.rows(), Observations.rows());
	
	covFunc.computeSymmetric(Sigma, Locations);

	mat Cpred;
	Cpred.set_size(Locations.rows(), Xpred.rows(), false);

	covFunc.computeCovariance(Cpred, Locations, Xpred);

	//mat alpha = C * Observations;

	covFunc.computeSymmetric(Sigma, Locations);
	
	cholSigma = computeCholesky(Sigma);

	vec alpha = ls_solve_chol(cholSigma, Observations);

	Mean = Cpred.transpose() * alpha;

	vec sigsq;

	sigsq.set_size(Xpred.rows(), false);
	
	covFunc.computeDiagonal(sigsq, Xpred);	

	Variance = sum(elem_mult((Cpred.transpose() * C), Cpred.transpose()), 2);	
	Variance = sigsq + Variance;	
	
}

double GaussianProcess::loglikelihood() const
{

	mat Sigma(Observations.size(), Observations.size());
	mat cholSigma(Observations.size(), Observations.size());

	covFunc.computeSymmetric(Sigma, Locations);

	cholSigma = computeCholesky(Sigma);

	mat invSigma = computeInverseFromCholesky(Sigma);
	vec alpha = invSigma * Observations;

//	vec alpha = ls_solve(Sigma, Observations);
	
	double out1 = 0.5 * dot(Observations,alpha);

	double out2 = sum(log(diag(cholSigma)));
		
	return out1 + out2 + 0.5*Observations.size()*log(2*pi);
	
}

vec GaussianProcess::getParametersVector() const
{
	return covFunc.getParameters();
}

void GaussianProcess::setParametersVector(const vec p)
{
	covFunc.setParameters(p);
}

double GaussianProcess::objective() const
{
	return loglikelihood();
}

/*
vec GaussianProcess::gradient() const
{
	vec grads(covFunc.getNumberParameters());

	
	mat Sigma(Observations.size(), Observations.size());
	mat cholSigma(Observations.size(), Observations.size());

	covFunc.computeSymmetric(Sigma, Locations);

	cholSigma = computeCholesky(Sigma);

	vec alpha = ls_solve(Sigma, Observations);
//	mat W = ls_solve(cholSigma, ls_solve(cholSigma.transpose(), eye(Observations.size()))) - outer_product(alpha, alpha, false);
	mat W = (inv(Sigma) - outer_product(alpha, alpha, false));

	mat partialDeriv(Observations.size(), Observations.size());

	for(int i = 0; i < covFunc.getNumberParameters(); i++)
	{
		covFunc.getParameterPartialDerivative(partialDeriv, i, Locations);
		//grads(i) = sum(sum(elem_mult(W, partialDeriv))) / 2;
		grads(i) = elem_mult_sum(W, partialDeriv) / 2;
// official - but slower		grads(i) = sum(diag(W * partialDeriv)) / 2;
	}
	return grads; 
}
*/

vec GaussianProcess::gradient() const
{
	vec grads(covFunc.getNumberParameters());

	
	mat Sigma(Observations.size(), Observations.size());
	mat cholSigma(Observations.size(), Observations.size());

	covFunc.computeSymmetric(Sigma, Locations);
	cholSigma = computeCholesky(Sigma);
	mat invSigma = computeInverseFromCholesky(Sigma);
	vec alpha = invSigma * Observations;




//	vec alpha = ls_solve(Sigma, Observations);
	mat W = (invSigma - outer_product(alpha, alpha, false));

	mat partialDeriv(Observations.size(), Observations.size());

	for(int i = 0; i < covFunc.getNumberParameters(); i++)
	{
		covFunc.getParameterPartialDerivative(partialDeriv, i, Locations);
		//grads(i) = sum(sum(elem_mult(W, partialDeriv))) / 2;
		grads(i) = elem_mult_sum(W, partialDeriv) / 2;
// official - but slower		grads(i) = sum(diag(W * partialDeriv)) / 2;
	}
	return grads; 
}





vec GaussianProcess::getGradientVector() const
{
	return gradient();
}

mat GaussianProcess::computeCholesky(const mat& iM) const 
{
	mat M = iM;// oops, was i inadvertantly writing to this?
	assert(M.rows() == M.cols());
	
	const double ampl = 1.0e-10;
	const int maxAttempts = 10;
	
	mat cholFactor(M.rows(), M.cols());

	int l = 0;
	bool success = chol(M, cholFactor);
	if(success)
	{
		return cholFactor;
	}
	else
	{
		double noiseFactor = abs(ampl * (trace(M) / double(M.rows())));
		while(!success)
		{
			M = M + (noiseFactor * eye(M.rows()));

			if(l > maxAttempts)
			{
				cerr << "Unable to compute cholesky decomposition" << endl;
				break;
			}
			l++;
			noiseFactor = noiseFactor * 10;
			success = chol(M, cholFactor);
		}
		cout << "Matrix not positive definite.  After " << l << " attempts, " << noiseFactor << " added to the diagonal" << endl;
	}
	return cholFactor;
}

mat GaussianProcess::computeInverseFromCholesky(const mat& C) const
{
	mat cholFactor = computeCholesky(C);
	mat invChol = backslash(cholFactor, eye(cholFactor.rows()));
	return invChol * invChol.transpose();
}

