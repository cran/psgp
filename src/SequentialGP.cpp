#include "SequentialGP.h"



SequentialGP::SequentialGP(int Inputs, int Outputs, mat& Xdata, vec& ydata, CovarianceFunction& cf) : Locations(Xdata), Observations(ydata), covFunc(cf), ForwardModel(Inputs, Outputs)
{
	assert(Locations.rows() == Observations.size());
	//covFunc = &cf;
	
	C.set_size(0, 0);
	KB.set_size(0, 0);
	Q.set_size(0, 0);
	Alpha.set_size(0, 1);
	ActiveSet.set_size(0, Inputs);
	
	sizeActiveSet = 0;
	
	// default value TODO
	maxActiveSet = 300;
	epsilonTolerance = 1e-5;
	
}

SequentialGP::~SequentialGP()
{
}

// train with one likelihood model
void SequentialGP::computePosterior(const LikelihoodType& noiseModel)
{
	assert(Locations.rows() == Observations.length());

	// randperm
	vec rndNums = itpp::randu(Locations.rows());
	ivec rndIdx = itpp::sort_index(rndNums);
	
	for(int i=0; i<Observations.length(); i++)	
	{
		int idx = rndIdx(i);
		//idx = i;
		// put: idx = i; for sequential index
		addOne(Locations.get_row(idx), Observations(idx), noiseModel);
	}
}

// train with multiple likelihood models
void SequentialGP::computePosterior(const ivec& LikelihoodModel, const Vec<LikelihoodType *> noiseModels)
{
	assert(Locations.rows() == LikelihoodModel.length());
	assert(Locations.rows() == Observations.length());

	// randperm
	vec rndNums = itpp::randu(Locations.rows());
	ivec rndIdx = itpp::sort_index(rndNums);
	
	for(int i=0; i<Observations.length(); i++)
	{
		int idx = rndIdx(i);
		//idx = i;
		// put: idx = i; for sequential index
		assert((LikelihoodModel(idx)) <= noiseModels.length());
		assert((LikelihoodModel(idx)) >= 1);
		addOne(Locations.get_row(idx), Observations(idx), *noiseModels(LikelihoodModel(idx) - 1));
	}
}


void SequentialGP::addOne(const vec& X, const double Observation, const LikelihoodType& noiseModel)
{
	mat KX, eHat, temp;
	mat Xmat(1, length(X));
	Xmat.set_row(0, X);	

	double sig0, sigx, rtp1, qtp1, mu, gamma;
	
	temp.set_size(1,1, false);
	covFunc.computeSymmetric(temp, Xmat);
	sig0 = temp(0,0);
	
	if(sizeActiveSet == 0)
	{
		sigx = sig0;
		mu = 0;
	}
	else
	{
		KX.set_size(sizeActiveSet, 1, false);
		covFunc.computeCovariance(KX, ActiveSet, Xmat);
		
		temp.set_size(sizeActiveSet,sizeActiveSet);
		temp = (KX.transpose() * C) * KX;
		temp = temp + sig0;
		sigx = temp(0,0);
		
		if(sigx < 1e-12)
		{
			sigx = 1e-12;
		//	cout << "sorting out variance" << endl;
		}
		
		temp.set_size(1, 1, false);
		temp = KX.transpose() * Alpha;
		mu = temp(0,0);		
	}
	
	noiseModel.updateCoefficients(qtp1, rtp1, Observation, mu, sigx);
	//cout << noiseModel.getID() << "( "<<qtp1<<", "<<rtp1<<")  ";
	//    [K1, K2]  = ogpparadj(K1,K2,cM+pM,cV,1e4./ikk,max(-1./ikk,-1e-3));

	//[tV,tu] = eig( - (eye(nout)+K2*sigX2)\K2 );
	//tu      = diag(tu);
	//iSmall  = find(tu<tolS);
	//iLarge  = find(tu>tolL);
	double tolL = 1e3;
	double tolS = 1e-6;
	double sqrtPt = sqrt(sigx);
	double tu = -sqrtPt*rtp1*sqrtPt; 
	bool mod = false;
	if(tu > tolL)
	{
		tu = tolL;
		cout << "Shrinking lambda" << endl;
		mod = true;
	}

	if(tu < tolS)
	{
		tu = tolS;
		cout << "(REV) Shrinking lambda" << endl;
		mod = true;
	}

	if(mod) {
		rtp1 = - (tu/sqrtPt)/tu;
		rtp1 = rtp1+itpp::eps;
		rtp1 = rtp1 + rtp1;
	}
	
//	cout << "   ( "<<qtp1<<", "<<rtp1<<") "<< Observations(0,0) <<", " << mu << ", "<< sigx << endl;
	
	if(sizeActiveSet == 0)
	{
		eHat = zeros(0, 1);
		gamma = sig0;
	}
	else
	{
		temp.set_size(sizeActiveSet, sizeActiveSet, false);
		eHat = inv(KB) * KX; // should use invKB
		//cout << "setting temp size" << endl;
		temp.set_size(1, 1, false);
		//cout << "KX'*eHat" << endl;
		temp = KX.transpose() * eHat;
		gamma = sig0 - temp(0,0);
	}
	
	if((gamma < epsilonTolerance) | (sizeActiveSet >= maxActiveSet))
	{	
		// Project
//		cout << "Sparse update" << endl;
		sparseUpdate(KX, eHat, qtp1, rtp1);
	}
	else
	{
		// add Active Point
//		cout << "Full update" << endl;
		fullUpdate(X, KX, eHat, gamma, qtp1, rtp1, sig0);
		// pruning TODO
	}

	
}

void SequentialGP::sparseUpdate(mat& KX, mat& eHat, const double qtp1, const double rtp1)
{
	mat sHat;
	
	// cout << "Sparse update" << endl;
	sHat = C * KX;
	sHat = sHat + eHat;
	Alpha = Alpha + (sHat * qtp1);
	C = C + ((sHat * sHat.transpose()) * rtp1);	
}

void SequentialGP::fullUpdate(const mat& X, const mat& KX, mat& eHat, const double gamma, const double qtp1, const double rtp1, const double sig0)
{
	mat stp1;

	 // cout << "Full update" << endl;
	
	sizeActiveSet = sizeActiveSet + 1;
	
	// Increase active set size
	ActiveSet.set_size(sizeActiveSet, Locations.cols(), true);

	ActiveSet.set_row(sizeActiveSet - 1, X.transpose().get_row(0));
			       	
	// update KB matrix
	KB.set_size(sizeActiveSet, sizeActiveSet, true);
	KB(sizeActiveSet - 1, sizeActiveSet - 1) = sig0;
	KB.set_submatrix(sizeActiveSet - 1, 0, KX.transpose());
	KB.set_submatrix(0, sizeActiveSet - 1, KX);

	// increase St+1 vector
	if(sizeActiveSet > 1)
	{
		stp1 = C * KX;
	}
	stp1.set_size(sizeActiveSet, 1, true);
	stp1(sizeActiveSet - 1, 0) = 1.0;

	// increase size of C and Alpha
	// then add new observation
	Alpha.set_size(sizeActiveSet, X.cols(), true);
	C.set_size(sizeActiveSet, sizeActiveSet, true);
	Alpha = Alpha + (stp1 * qtp1);
	C = C + ((stp1 * stp1.transpose()) * rtp1);
	
	// update Q matrix
	eHat.set_size(sizeActiveSet, 1, true);
	eHat(sizeActiveSet - 1, 0) = -1.0;
	Q.set_size(sizeActiveSet, sizeActiveSet, true);
	Q = Q + ((1.0 / gamma) * (eHat * eHat.transpose()));	
}


void SequentialGP::makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const
{
	assert(Mean.length() == Variance.length());
	assert(Xpred.rows() == Mean.length());
	
	// mean predictions
	mat Cpred(Xpred.rows(), ActiveSet.rows());
//	Cpred.set_size(Xpred.rows(), ActiveSet.rows(), false);
	covFunc.computeCovariance(Cpred, Xpred, ActiveSet);
	Mean = Cpred * Alpha;
	// variance predictions
	vec sigsq(Xpred.rows());
//	sigsq.set_size(Xpred.rows(), false);
	covFunc.computeDiagonal(sigsq, Xpred);
//	sigsq = itpp::diag(sigsq);	
	
//	Variance.set_col(0, itpp::sum(itpp::elem_mult((Cpred * C), Cpred), 2));
//	Variance = sigsq + Variance;
	Variance = sigsq + itpp::sum(itpp::elem_mult((Cpred * C), Cpred), 2);
}

vec SequentialGP::getParametersVector() const
{
	vec a = "0.0";
	return a;
}

void SequentialGP::setParametersVector(const vec p)
{

}

double SequentialGP::objective() const
{
	mat Z = KB * C;

	mat Sigma(Observations.size(), Observations.size());
	mat cholSigma(Observations.size(), Observations.size());

	cholSigma = computeCholesky(Z);

	mat invSigma = computeInverseFromCholesky(cholSigma);

	vec alpha = Alpha.get_col(0);

//	vec lP = alpha * invSigma;
	vec lP = alpha;
	vec rP = KB * alpha;

	double out1 = dot(lP , rP);
	double out2 = sum(log(diag(cholSigma)));
	
	return out1 + out2 + 0.5*Observations.size()*log(2*pi);

}

vec SequentialGP::gradient() const
{
	// TODO: implement this function
	vec a = "0.0";
	
	return a;
}

void SequentialGP::estimateParameters()
{
	// TODO: implement this function	
}

void SequentialGP::displayCovarianceParameters() const
{
	cout << "Model: Sequential Gaussian Process" << endl;
	cout << "  Active set size: " << ActiveSet.rows() << endl;
	cout << "  Kernel Matrix size: " << KB.rows() << " x " << KB.cols() << endl;
	cout << "  Inverse Kernel Matrix size: " << Q.rows() << " x " << Q.cols() << endl;
	cout << "  Alpha size: " << Alpha.rows() << " x " << Alpha.cols() << endl;
	cout << "  C size: " << C.rows() << " x " << C.cols() << endl;

}

mat SequentialGP::computeCholesky(const mat& iM) const 
{
	mat M = iM;
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

mat SequentialGP::computeInverseFromCholesky(const mat& C) const
{
	mat cholFactor = computeCholesky(C);
	mat invChol = backslash(cholFactor, eye(cholFactor.rows()));
	return invChol * invChol.transpose();
}

