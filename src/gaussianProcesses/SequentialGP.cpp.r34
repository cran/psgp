#include "SequentialGP.h"

////////////////////////////////////////////////////////////////////////////////
//
// CONSTRUCTOR
//
////////////////////////////////////////////////////////////////////////////////
SequentialGP::SequentialGP(int Inputs, int Outputs, int nActivePoints, mat& Xdata, vec& ydata, 
		                   CovarianceFunction& cf, int _iterChanging) : ForwardModel(Inputs, Outputs), 
		                   Locations(Xdata), Observations(ydata), covFunc(cf)
{
	assert(Locations.rows() == Observations.size());
	
	C.set_size(0, 0);
	KB.set_size(0, 0);
	Q.set_size(0, 0);
	Alpha.set_size(0);
	ActiveSet.set_size(0, Inputs);
	idxActiveSet.set_size(0);

	projectionP = zeros(Observations.length(), 0);
	lambdaP = zeros(Observations.length());
	alphaP = zeros(Observations.length());
	logZ = zeros(Observations.length());
	
	sizeActiveSet = 0;
	
	// default value TODO
	maxActiveSet = nActivePoints;
	epsilonTolerance = 1e-6;

	momentProjection = true;
	iterChanging = _iterChanging;
	iterFixed = 1;

//	likelihoodType = FullEvid;
	likelihoodType = UpperBound;
//	likelihoodType = Approximate;
	
	// If this is set to true, the locations already in the active set
	// will not be added+deleted again.
	selectiveSweep = false;
}

////////////////////////////////////////////////////////////////////////////////
//
// DESTRUCTOR
//
////////////////////////////////////////////////////////////////////////////////
SequentialGP::~SequentialGP()
{

}

// Toggle the "selective sweep" option
void SequentialGP::setSelectiveSweep(bool b) { selectiveSweep = b; }

////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE POSTERIOR WITH SINGLE LIKELIHOOD MODEL
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::computePosterior(const LikelihoodType& noiseModel)
{
	assert(Locations.rows() == Observations.length());

	bool fixActiveSet = false;


	for(int cycle = 1; cycle <= (iterChanging + iterFixed); cycle++)
	{
		if(cycle > iterChanging)
		{
			fixActiveSet = true;
		}

		// do a randperm
		vec rndNums = randu(Locations.rows());
		ivec rndIdx = sort_index(rndNums);

		// Only add locations which are not already in the active set
		/*
		bvec activeIndex(Observations.length());
		
		if (selectiveSweep) {
    		for (int i=0; i<activeIndex.length(); i++)
    		    activeIndex(i) = true;
    		for (int i=0; i<idxActiveSet.length(); i++) 
    		    activeIndex(idxActiveSet(i)) = false;
		}
		*/
		
		// send the observations individually
		
		int nobs = Observations.length();
		for(int i=0; i<nobs; i++)	
		{
		    // If selectiveSweep is set to true, only add location if it is not already 
		    // in the active set. Otherwise, always add observation (even if already in
		    // active set).
		    /*
		    if (selectiveSweep) {   
		        if (activeIndex(rndIdx(i))) 
		            addOne(rndIdx(i), noiseModel, fixActiveSet);
		        else
		            cout << "Skipping observation " << rndIdx(i) << " (already in active set)" << endl; 
		    }
		    else
		    */ 
		    
		    // cout << "\rAdding observations: " << rndIdx(i) << "/" << nobs;
		    addOne(rndIdx(i), noiseModel, fixActiveSet);
		}
		// cout << endl;

	}

}



////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE POSTERIOR WITH SINGLE LIKELIHOOD MODEL
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::computePosteriorFixedActiveSet(const LikelihoodType& noiseModel, ivec iActive)
{
	assert(Locations.rows() == Observations.length());
	
	bool fixActiveSet = false;

	// Initialise active set
	for(int i=0; i<iActive.length(); i++)  
	{
	    cout << "\rAdding observations: " << i << "/" << iActive.length();
	    addOne(iActive(i), noiseModel, fixActiveSet);
	}
	cout << endl;
	recomputePosteriorFixedActiveSet(noiseModel);
}

/**
 * Revisit data with fixed active set
 */
void SequentialGP::recomputePosteriorFixedActiveSet(const LikelihoodType& noiseModel)
{
	// Do a randperm
	vec rndNums = randu(Locations.rows());
	ivec rndIdx = sort_index(rndNums);

	// Visit all observations
	bool fixActiveSet = true;
	for(int i=0; i<Observations.length(); i++)	
	{
	    cout << "\rAdding observations: " << i << "/" << Observations.length();
	    addOne(rndIdx(i), noiseModel, fixActiveSet);
	}
	cout << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE POSTERIOR WITH MULTIPLE LIKELIHOODS
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::computePosterior(const ivec& LikelihoodModel, const Vec<LikelihoodType *> noiseModels)
{
	assert(Locations.rows() == LikelihoodModel.length());
	assert(Locations.rows() == Observations.length());

	bool fixActiveSet = false;
	
	for(int cycle = 1; cycle <= (iterChanging + iterFixed); cycle++)
	{
		if(cycle > iterChanging)
		{
			fixActiveSet = true;
		}

		// randperm
		vec rndNums = randu(Locations.rows());
		ivec rndIdx = sort_index(rndNums);

		for(int i=0; i<Observations.length(); i++)	
		{
			assert((LikelihoodModel(rndIdx(i))) <= noiseModels.length());
			assert((LikelihoodModel(rndIdx(i))) >= 1);
			addOne(rndIdx(i), *noiseModels(LikelihoodModel(rndIdx(i)) - 1), fixActiveSet);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// ADD ONE OBSERVATION
//
////////////////////////////////////////////////////////////////////////////////
inline void SequentialGP::addOne(int index, const LikelihoodType& noiseModel, const bool fixActiveSet)
{
	mat KX, temp;
	mat Xmat = Locations.get_row(index).transpose();

	double sig0, sigx, rtp1, qtp1, mu, gamma, logEvidence;
	double Observation = Observations(index);
	vec eHat;	

//	cout << index << ", ";

	temp.set_size(1,1, false);
	covFunc.computeSymmetric(temp, Xmat);

	sig0 = temp(0,0);
	
	// SITE REMOVAL
	if(abs(lambdaP(index)) > 1e-10)
	{
		// h_i = Sigma*Phi_i
		vec projectionRow = projectionP.get_row(index);
		vec teK = KB * projectionRow;
		
		// Eq. 4.19 and 4.20 in Lehel
		vec tMe = (eye(sizeActiveSet) + (C * KB)) * projectionRow;
		double tV = lambdaP(index) / (1 - (lambdaP(index) * dot(teK, tMe))); // do we need a check for a zero here?
		Alpha = Alpha + tMe * tV * (dot(teK, Alpha) - alphaP(index));
		C = C + (outer_product(tMe, tMe) * tV);
	}

	// CALCULATE CAVITY MEAN/VARIANCE
	if(sizeActiveSet == 0)
	{
		sigx = sig0;
		mu = 0;
	}
	else
	{
		KX.set_size(sizeActiveSet, 1, false);
		covFunc.computeCovariance(KX, ActiveSet, Xmat);
		
		// sigx = predictive variance 
		temp = (KX.transpose() * C) * KX;
		sigx = sig0 + temp(0,0);
		
		// mu = predictive mean
		temp = KX.transpose() * Alpha;
		mu = temp(0,0);
	}

	// STABILITY
	if(sigx < 1e-12)
	{
		sigx = 1e-12;
	}

//cout << "Noise model" << endl;
	// UPDATE COEFFICIENTS FOR SPECIFIED LIKELIHOOD
	// Lehel: Sec. 2.4, Eq. 2.42, 3.28
	logEvidence = noiseModel.updateCoefficients(qtp1, rtp1, Observation, mu, sigx);


	// stablize update coefficients
	stabiliseCoefficients(qtp1, rtp1, mu, sigx, 1e3, 1e-6);


////////////////////////////////////////////////////////////////////////////////
	
	if(sizeActiveSet == 0)
	{
		eHat = zeros(0);
		gamma = sig0;
	}
	else
	{
		temp.set_size(sizeActiveSet, sizeActiveSet, false);

////////// LOOK FOR MORE EFFICIENT ALTERNATIVES
//		eHat = inv(KB) * KX; // should use invKB
		
		// Lehel: Eq. 3.3
		eHat = backslash(KB, KX);

		temp = KX.transpose() * eHat;
		gamma = sig0 - temp(0,0);
	}
	
	// STABILITY
	if(gamma < 1.0e-12)
	{
		gamma = 0.0;
	}


	
	// FULL OR SPARSE MODEL UPDATE
	if((gamma < epsilonTolerance) | fixActiveSet)
	{	
//        cout << "Update sparse" << endl;
		updateSparse(KX, eHat, gamma, qtp1, rtp1, mu, sigx, logEvidence, index);
	}
	else
	{
//		cout << "Update full" << endl;
		updateFull(KX, eHat, gamma, qtp1, rtp1, sig0, mu, sigx, logEvidence, index);

//cout << "Extra" << endl;
		// REMOVE EXTRA ACTIVE POINTS
		if(!fixActiveSet)
		{
			while(sizeActiveSet > maxActiveSet)
			{
				int removalCandidate = findLeastInformativeActivePoint(FullKL);
				// cout << "Deleting active point " << idxActiveSet(removalCandidate) << endl;
				deleteActivePoint(removalCandidate);
			}			
		}

//cout << "Geometry" << endl;
		// REMOVE ACTIVE POINTS BASED ON GEOMETRY
		// For stability, if variance == 0 remove associated active point 
		while(sizeActiveSet > 0)
		{
			int removalCandidate = findLeastInformativeActivePoint(Geometric);
			if(removalCandidate == -1)
			{
				break;
			}
			else
			{
				cout << "Removing point (covariance collapsed)" << endl;
				deleteActivePoint(removalCandidate);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// SPARSE UPDATE
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::updateSparse(mat& KX, vec& eHat, const double gamma, const double qtp1, 
		                        const double rtp1, const double currentMean, const double currentVar, 
		                        const double logEvidence, const int index)
{
	double ratio;
	// Lehel: Eq. 3.4
	vec sHat = (C * KX).get_col(0) + eHat;
	mat sHatOuter = outer_product(sHat, sHat, false);

	// Lehel: Eq. 3.29
	// Check if momentProjection needs to be toggled??
	if(momentProjection)
	{
		double eta = 1.0 + (rtp1 * gamma);
		Alpha = Alpha + (sHat * (qtp1 / eta));
		C = C + (sHatOuter * (rtp1 / eta));
	}
	else
	{
		// Lehel: Sec. 3.5 
		Alpha = Alpha + (sHat * qtp1);
		C = C + (sHatOuter * rtp1);
	}

	assert(eHat.length() == projectionP.cols());
	assert(sizeActiveSet == projectionP.cols());

	// update EP parameters
	ratio = qtp1 / rtp1;
	logZ(index)= logEvidence + (log(2 * pi) - log(abs(rtp1)) - (qtp1 * ratio)) / 2.0;
	projectionP.set_row(index, eHat);
	alphaP(index) = currentMean - ratio;
	lambdaP(index) = rtp1 / (1.0 + (rtp1 * currentVar));

	if(!momentProjection)
	{
		lambdaP(index) = lambdaP(index) / (1.0 + (gamma * lambdaP(index)));
	}

}


	
////////////////////////////////////////////////////////////////////////////////
//
// FULL UPDATE
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::updateFull(const mat& KX, vec& eHat, const double gamma, const double qtp1, const double rtp1, const double sig0, const double currentMean, const double currentVar, const double logEvidence, const int index)
{
	double ratio;
	vec stp1;

	sizeActiveSet = sizeActiveSet + 1;
	
	idxActiveSet.set_size(sizeActiveSet, true);
	idxActiveSet(sizeActiveSet - 1) = index;

	// Increase active set size
	ActiveSet.set_size(sizeActiveSet, Locations.cols(), true);
	ActiveSet.set_row(sizeActiveSet - 1, Locations.get_row(index));

	// update KB matrix
	KB.set_size(sizeActiveSet, sizeActiveSet, true);
	KB(sizeActiveSet - 1, sizeActiveSet - 1) = sig0;
	KB.set_submatrix(sizeActiveSet - 1, 0, KX.transpose());
	KB.set_submatrix(0, sizeActiveSet - 1, KX);

	// increase St+1 vector
	if(sizeActiveSet > 1)
	{
		stp1 = (C * KX).get_col(0);
	}
	stp1.set_size(sizeActiveSet, true);
	stp1(sizeActiveSet - 1) = 1.0;

	// increase size of C and Alpha
	// then add new observation
	Alpha.set_size(sizeActiveSet, true);
	C.set_size(sizeActiveSet, sizeActiveSet, true);
	Alpha = Alpha + (stp1 * qtp1);
	C = C + (outer_product(stp1, stp1, false) * rtp1);	

	// update Q matrix
	eHat.set_size(sizeActiveSet, true);
	eHat(sizeActiveSet - 1) = -1.0;
	Q.set_size(sizeActiveSet, sizeActiveSet, true);
	Q = Q + (outer_product(eHat, eHat, false) / gamma);	

	// update EP parameters
	ratio = qtp1 / rtp1;
	logZ(index)= logEvidence + (log(2 * pi) - log(abs(rtp1)) - (qtp1 * ratio)) / 2.0;
	
	// Double check here (Lehel: Eq. 4.30)
	projectionP.set_size(Observations.length(), sizeActiveSet, true);
	projectionP.set_row(index, zeros(sizeActiveSet));
	projectionP(index, sizeActiveSet - 1) = 1;

	alphaP(index) = currentMean - ratio;
	lambdaP(index) = rtp1 / (1.0 + (rtp1 * currentVar));

}

////////////////////////////////////////////////////////////////////////////////
//
// SCORE ACTIVE POINTS
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::scoreActivePoints(ScoringMethod sm)
{
	vec diagC, diagS, term1, term2, term3;
	vec diagInvGram = diag(Q);


	switch(sm)
	{
		case Geometric : return(1.0 / diagInvGram);
			break;
		case MeanComponent : return(elem_div(elem_mult(Alpha,Alpha), diag(C) + diagInvGram));
			break;
		case FullKL : // Lehel: Eq. 3.23 
			diagC = diag(C);
			diagS = diag((projectionP.transpose() * diag(lambdaP)) * projectionP); // need to find a more efficient way of doing this
			term1 = elem_div(elem_mult(Alpha,Alpha), diagC + diagInvGram);
			term2 = elem_div(diagS, diagInvGram);
			term3 = log(1.0 + elem_div(diagC , diagInvGram));
			return (term1 + term2 + term3);
			break;
		default : cerr << "Unknown scoring method" << endl;
			break;
	}
	return zeros(sizeActiveSet);
}

////////////////////////////////////////////////////////////////////////////////
//
// FIND LEAST INFORMATIVE ACTIVE POINT
//
////////////////////////////////////////////////////////////////////////////////
int SequentialGP::findLeastInformativeActivePoint(ScoringMethod sm)
{
	vec scores = scoreActivePoints(sm);
	int leastInformativeIndex = min_index(scores);

	switch(sm)
	{
		case Geometric :
			// check to see if the informativeness is above a certain point
			if(scores(leastInformativeIndex) > (epsilonTolerance / 1000.0))
			{
				return(-1);
			}
			break;

		case MeanComponent :
		case FullKL :
		default :
			// do nothing
			break;
	}

	return leastInformativeIndex;
}

////////////////////////////////////////////////////////////////////////////////
//
// DELETE ACTIVE POINT
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::deleteActivePoint(int index)
{
	assert(sizeActiveSet == Q.rows());
	double q_star, c_star;
	vec red_q(sizeActiveSet - 1), red_Ag(sizeActiveSet - 1);

	// calculate set difference
	int delIdx = index;
	ivec restIdx(sizeActiveSet - 1);
	int temp = 0;
	for(int i = 0; i < sizeActiveSet; i++)
	{
		if(i != index)
		{
			restIdx(temp++) = i;
		}
	}	

	// this bit could be considered to be a bit messy
	q_star = Q(delIdx, delIdx);
	c_star = C(delIdx, delIdx);

	// red_q  = net.KBinv(exI,rmI);
	for(int i = 0; i < (sizeActiveSet - 1); i++)
	{
		red_q(i) = Q(restIdx(i), delIdx);			
	}
	
	if(momentProjection)
	{
		// These correspond to Lehel Appendix D
		// Unsure what the other case corresponds to
		
		// red_Ag = red_q + net.C(exI,rmI);
		for(int i = 0; i < (sizeActiveSet - 1); i++)
		{
			red_Ag(i) = red_q(i) + C(restIdx(i), delIdx);			
		}
		double delAlpha = Alpha(delIdx);
		Alpha.del(delIdx);
		Alpha = Alpha - (red_Ag * (delAlpha / (q_star + c_star)));

		C.del_row(delIdx);
		C.del_col(delIdx);
		C = C + outer_product(red_q, red_q / q_star) - outer_product(red_Ag, red_Ag / (q_star + c_star));		

		//  net.w  = net.w(exI,:) - red_Ag * ((q_star+c_star)\net.w(rmI));
		//  net.C  = net.C(exI,exI) + red_q * (q_star\red_q') - red_Ag * ((q_star+c_star)\red_Ag');

	}
	else
	{
		
		vec tempQ = red_q / q_star; //  tempQ = red_q/q_star;
		mat matQ;
		double delAlpha = Alpha(delIdx);
		Alpha.del(delIdx);
		Alpha = Alpha - (tempQ * delAlpha); //  net.w = net.w(exI) - tempQ*net.w(rmI);

		// red_c = net.C(rmI,exI);
		vec red_c(sizeActiveSet - 1);
		for(int i = 0; i < (sizeActiveSet - 1); i++)
		{
			red_c(i) = C(delIdx, restIdx(i));			
		}
		C.del_row(delIdx);
		C.del_col(delIdx);

		C = C + outer_product(tempQ * c_star, tempQ); //  net.C = net.C(exI,exI) + tempQ*c_star*tempQ';
		matQ = outer_product(tempQ, red_c); //  tempQ = tempQ*red_c;
		C = (C - matQ) - matQ.transpose();	//  net.C = net.C - tempQ - tempQ';
	}


	if(!momentProjection)
	{	
		vec prDel = pow(projectionP.get_col(delIdx), 2.0);
		vec lDel = elem_mult(lambdaP, prDel);
		lambdaP = lambdaP + (elem_div(pow(lDel, 2.0), q_star + lDel));	

		// prDel = ep.projP(:,iBV).^2;
		// lDel  = full(ep.lamP*prDel);
		// ep.lamP = ep.lamP + diag(lDel.^2./(q_star + lDel));
  	}
	
	// ep.projP = ep.projP(:,exI) - ep.projP(:,rmI) * (q_star \ red_q');
	vec projDel = projectionP.get_col(delIdx);
	projectionP.del_col(delIdx);
	projectionP = projectionP - outer_product(projDel, red_q / q_star);


	Q.del_row(delIdx);
	Q.del_col(delIdx); // net.KBinv  = net.KBinv(exI,exI) - red_q*(q_star\red_q');
	Q = Q - (outer_product(red_q, red_q) / q_star);
	

	KB.del_row(delIdx);
	KB.del_col(delIdx); // net.KB   = net.KB(exI,exI);
	ActiveSet.del_row(delIdx); // net.BV   = net.BV(exS,:);

	idxActiveSet.del(delIdx);
	sizeActiveSet--;
}

////////////////////////////////////////////////////////////////////////////////
//
// CHECK UPDATE COEFFICIENTS TO ENSURE STABILITY
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::stabiliseCoefficients(double& q, double& r, double cavityMean, double cavityVar, double upperTolerance, double lowerTolerance)
{

	double sqrtPt = sqrt(cavityVar);
	double tu = -sqrtPt * r * sqrtPt; 
	bool mod = false;
	if(tu > upperTolerance)
	{
		tu = upperTolerance;
		cout << "Shrinking lambda" << endl;
		mod = true;
	}

	if(tu < lowerTolerance)
	{
		tu = lowerTolerance;
		cout << "(REV) Shrinking lambda" << endl;
		mod = true;
	}

	if(mod) {
		r = - (tu / sqrtPt) / tu;
		r = r + itpp::eps;
		r = r + r;
	}
}


////////////////////////////////////////////////////////////////////////////////
//
// MAKE PREDICTIONS
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::makePredictions(vec& Mean, vec& Variance, const mat& Xpred, 
		                           CovarianceFunction& cf) const
{
	assert(Mean.length() == Variance.length());
	assert(Xpred.rows() == Mean.length());
	
	mat Cpred(Xpred.rows(), ActiveSet.rows());
	cf.computeCovariance(Cpred, Xpred, ActiveSet);

	Mean = Cpred * Alpha;
	vec sigsq(Xpred.rows());
	covFunc.computeDiagonal(sigsq, Xpred);
	Variance = sigsq + sum(elem_mult((Cpred * C), Cpred), 2);
}

/**
 * Same as above, but using the current (stored) covariance function to make
 * the predictions.
 **/
void SequentialGP::makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const
{
	makePredictions(Mean, Variance, Xpred, covFunc);
}

////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE SIMULATIONS
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::simulate(const mat& Xpred, bool approx) const
{
	mat cov, vCov, kxbv(Xpred.rows(), sizeActiveSet);
	vec dCov, samp;

	covFunc.computeCovariance(kxbv, Xpred, ActiveSet);

	if(approx)
	{
		cov = Q + C;
		vCov.set_size(sizeActiveSet, sizeActiveSet);
		dCov.set_size(sizeActiveSet);
		eig_sym(cov, dCov, vCov);
		vCov = kxbv * vCov;
		samp = randn(sizeActiveSet);
	}
	else
	{
		mat kxx(Xpred.rows(), Xpred.rows());
		covFunc.computeSymmetric(kxx, Xpred);
		cov = kxx + ((kxbv * C) * kxbv.transpose());
		eig_sym(cov, dCov, vCov);
		samp = randn(Xpred.rows());
	}

	dCov = sqrt(abs(dCov));

	vec a1 = kxbv * Alpha;
	mat a2 = vCov * diag(dCov);

	return(a1 + (a2 * samp));
}


////////////////////////////////////////////////////////////////////////////////
//
// GET MODEL PARAMETERS
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::getParametersVector() const
{
	return covFunc.getParameters();
}

////////////////////////////////////////////////////////////////////////////////
//
// SET MODEL PARAMETERS
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::setParametersVector(const vec p)
{
//	cout << "Set parameters vector"<< endl;
	covFunc.setParameters(p);
}

////////////////////////////////////////////////////////////////////////////////
//
// RECOMPUTE POSTERIOR PARAMETERS
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::recomputePosterior()
{
	mat KBold = KB;
	mat Kplus(Observations.length(), sizeActiveSet);

	covFunc.computeSymmetric(KB, ActiveSet);
	covFunc.computeCovariance(Kplus, Locations, ActiveSet);
	backslash(KB.transpose(), Kplus.transpose(), projectionP); // output is rightmost parameter
	mat lambdaPdiag = diag(lambdaP);
	mat projLam = projectionP * lambdaPdiag;
	mat UU = projLam * projectionP.transpose();
	mat CC = UU * KB + eye(sizeActiveSet);
	Alpha = backslash(CC, projLam * alphaP);
	C = -backslash(CC, UU);
	Q = computeInverseFromCholesky(KB);
}

void SequentialGP::updateModel() {
    recomputePosterior();
}

////////////////////////////////////////////////////////////////////////////////
//
// RESET POSTERIOR REPRESENTATION
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::resetPosterior()
{
	C.set_size(0, 0);
	KB.set_size(0, 0);
	Q.set_size(0, 0);
	Alpha.set_size(0);
	ActiveSet.set_size(0, getInputDimensions());
	idxActiveSet.set_size(0);
	projectionP = zeros(Observations.length(), 0);
	lambdaP = zeros(Observations.length());
	alphaP = zeros(Observations.length());
	logZ = zeros(Observations.length());

	sizeActiveSetOld = sizeActiveSet;
	sizeActiveSet = 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// OBJECTIVE FUNCTION
//
////////////////////////////////////////////////////////////////////////////////
double SequentialGP::objective() const
{
	double evidence;

	switch(likelihoodType)
	{
		case FullEvid : evidence = compEvidence();
			break;
		case Approximate : evidence = compEvidenceApproximate();
			break;
		case UpperBound :	evidence = compEvidenceUpperBound();
			break;
		default:	return 0.0;
	}

//	cout << "Evidence: " << evidence << endl;

	return evidence;
}

////////////////////////////////////////////////////////////////////////////////
//
// GRADIENT FUNCTION
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::gradient() const
{
	vec g;

	switch(likelihoodType)
	{
		case FullEvid : g = gradientEvidence();
			break;
		case Approximate : g = gradientEvidenceApproximate();
			break;
		case UpperBound :	g = gradientEvidenceUpperBound();
			break;
		default : g.set_size(covFunc.getNumberParameters());
	}

	return g;
}

////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE MODEL EVIDENCE
//
////////////////////////////////////////////////////////////////////////////////
double SequentialGP::compEvidence() const
{

	cvec es;
	mat KB_new(sizeActiveSet, sizeActiveSet);

	covFunc.computeSymmetric(KB_new, ActiveSet);

	double evid = sum(log(lambdaP));

	evid -= sum(elem_mult(pow(alphaP, 2.0), lambdaP));
	evid += 2.0 * sum(logZ);
	evid -= lambdaP.length() * log(2.0 * pi);
	mat Klp = projectionP.transpose() * diag(lambdaP);
	mat Ksm = (Klp * projectionP) * KB_new + eye(sizeActiveSet);
	vec Kall = Klp * alphaP;	
	mat Kinv= backslash(Ksm.transpose(), KB_new.transpose());

	evid += dot(Kall, Kinv.transpose() * Kall);

	if(eig(Ksm, es))
	{
		evid -= sum(log(real(es)));
	}
	else
	{
		cerr << "Problem computing evidence: eig_sum()" << endl;
	}

	return -evid / 2.0;
}

////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE APPROXIMATE MODEL EVIDENCE
//
////////////////////////////////////////////////////////////////////////////////
double SequentialGP::compEvidenceApproximate() const
{
	mat cholSigma(sizeActiveSet, sizeActiveSet);
	mat Sigma(sizeActiveSet, sizeActiveSet);
	covFunc.computeSymmetric(Sigma, ActiveSet);
	mat invSigma = computeInverseFromCholesky(Sigma);
	vec obsActiveSet = Observations(idxActiveSet);
	vec alpha = invSigma * obsActiveSet;

	double like1 = sum(log(diag(computeCholesky(Sigma))));
	double like2 = 0.5 * dot(obsActiveSet, alpha);

	return like1 + like2 + 0.5 * ActiveSet.size() * log(2 * pi);
}

////////////////////////////////////////////////////////////////////////////////
//
// COMPUTE EVIDENCE UPPER BOUND
//
////////////////////////////////////////////////////////////////////////////////
double SequentialGP::compEvidenceUpperBound() const
{
	mat KB_new(sizeActiveSet, sizeActiveSet);
	covFunc.computeSymmetric(KB_new, ActiveSet);
	/*
	double cond = itppext::cond(KB_new);
	while (cond < 1e-6) {
	    KB_new = KB_new + 1e-6*eye(KB_new.rows());
	    cout << "Ill-conditionned matrix (cond = " << cond << ")" << endl;
	    cond = itppext::cond(KB_new);
	}
	*/
	double like1 = 2.0 * (sum(log(diag(computeCholesky(KB_new)))));
	double like2 = trace((eye(sizeActiveSet) + 
	       (KB * (C + outer_product(Alpha, Alpha)))) * backslash(KB_new, KB));
	// cout << "lik1 = " << like1 << endl;
	// cout << "lik2 = " << like2 << endl;
	return like1 + like2;
}

////////////////////////////////////////////////////////////////////////////////
//
// GRADIENTS FOR EVIDENCE
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::gradientEvidence() const
{
	vec grads = zeros(covFunc.getNumberParameters());
	return grads;

}

////////////////////////////////////////////////////////////////////////////////
//
// GRADIENTS FOR APPROXIMATE EVIDENCE
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::gradientEvidenceApproximate() const
{
	vec grads(covFunc.getNumberParameters());

	mat cholSigma(sizeActiveSet, sizeActiveSet);
	mat Sigma(sizeActiveSet, sizeActiveSet);

	covFunc.computeSymmetric(Sigma, ActiveSet);
	cholSigma = computeCholesky(Sigma);
	mat invSigma = computeInverseFromCholesky(Sigma);
	vec obsActiveSet = Observations(idxActiveSet);
	vec alpha = invSigma * obsActiveSet;


	mat W = (invSigma - outer_product(alpha, alpha, false));

	mat partialDeriv(sizeActiveSet, sizeActiveSet);

	for(int i = 0; i < covFunc.getNumberParameters(); i++)
	{
		covFunc.getParameterPartialDerivative(partialDeriv, i, ActiveSet);
		grads(i) = elem_mult_sum(W, partialDeriv) / 2.0;
	}
	return grads; 
}

////////////////////////////////////////////////////////////////////////////////
//
// GRADIENT FOR UPPER BOUND ON EVIDENCE
//
////////////////////////////////////////////////////////////////////////////////
vec SequentialGP::gradientEvidenceUpperBound() const
{

	vec grads(covFunc.getNumberParameters());

	mat W = eye(sizeActiveSet);
	mat KB_new(sizeActiveSet, sizeActiveSet);
	covFunc.computeSymmetric(KB_new, ActiveSet);
	/*
	double cond = itppext::cond(KB_new);
	while (cond < 1e-6) {
	    KB_new = KB_new + 1e-6*eye(KB_new.rows());
	    cout << "Ill-conditionned matrix (cond = " << cond << ")" << endl;
	    cond = itppext::cond(KB_new);
	}
	*/
	// prob with lengthscale calculation somewhere here
	// W = backslash(KB_new, W - (W + (KB * (C + outer_product(Alpha, Alpha)))) * backslash(KB_new, KB));
	
	// This gives the correct gradient for the length scale
	mat partialDeriv(sizeActiveSet, sizeActiveSet);
	mat U = backslash(KB_new,KB);
	
	W =  W + (KB * (C + outer_product(Alpha, Alpha)));
	
	for(int i = 0; i < covFunc.getNumberParameters(); i++)
	{
		covFunc.getParameterPartialDerivative(partialDeriv, i, ActiveSet);
		mat V1 = backslash(KB_new,partialDeriv);
		mat V2 = W*backslash(KB_new,partialDeriv*U);
		 
		grads(i) = trace(V1-V2);
	}

	return grads;
}

void SequentialGP::setLikelihoodType(LikelihoodCalculation lc)
{
	likelihoodType = lc;
}



////////////////////////////////////////////////////////////////////////////////
//
// DISPLAY MODEL PARAMETERS
//
////////////////////////////////////////////////////////////////////////////////
void SequentialGP::displayModelParameters() const
{
	cout << "Summary Sequential Gaussian Process" << endl;
	cout << "  Kernel Matrix size         : " << KB.rows() << " x " << KB.cols() << endl;
	cout << "  Inverse Kernel Matrix size : " << Q.rows() << " x " << Q.cols() << endl;
	cout << "  Alpha size                 : " << Alpha.size() << endl;
	cout << "  C size                     : " << C.rows() << " x " << C.cols() << endl;
	cout << "  Projection matrix size     : " << projectionP.rows() << " x " << projectionP.cols() << endl;
	cout << "  Lambda                     : " << lambdaP.size() << endl;
	cout << "  projection alpha           : " << alphaP.size() << endl;
	cout << "  log evidence vector        : " << logZ.size() << endl;
	cout << "  ----------------------------" << endl;
	cout << "  Predicion locations        : " << Locations.rows() << " x " << Locations.cols() << endl;
	cout << "  Observations               : " << Observations.size() << endl;
	cout << "  Active set size            : " << ActiveSet.rows() << " (max = " << maxActiveSet << ")" << endl;
	cout << "  Epsilon tolerance          : " << epsilonTolerance << endl;
	cout << "  Iterations Changing/Fixed  : " << iterChanging << "/" << iterFixed << endl;
	cout << "  Moment projection          : " << momentProjection << endl;
}

////////////////////////////////////////////////////////////////////////////////
//
// Auxilary functions
//
////////////////////////////////////////////////////////////////////////////////

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

