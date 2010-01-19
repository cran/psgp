/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008-2009                                   *
 *                                                                         *
 *   Ben Ingram, IngramBR@Aston.ac.uk                                      *
 *   Neural Computing Research Group,                                      *
 *   Aston University,                                                     *
 *   Aston Street, Aston Triangle,                                         *
 *   Birmingham. B4 7ET.                                                   *
 *   United Kingdom                                                        *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef SEQUENTIALGP_H_
#define SEQUENTIALGP_H_

#include <itpp/itbase.h>

#include "ForwardModel.h"
#include "../optimisation/Optimisable.h"
#include "../covarianceFunctions/CovarianceFunction.h"
#include "../likelihoodModels/LikelihoodType.h"
#include "../likelihoodModels/GaussianSampLikelihood.h"
// #include "itppext/itppext.h"

#include <cassert>

using namespace std;
using namespace itpp;


#ifndef PSGP_H_
  enum ScoringMethod {Geometric, MeanComponent, FullKL};
  enum LikelihoodCalculation {FullEvid, Approximate, UpperBound};
#endif 

class SequentialGP : public ForwardModel, public Optimisable
{
public:
	SequentialGP(int Inputs, int Outputs, int nActivePoints, mat& Xdata, vec& ydata, CovarianceFunction& cf, int _iterChanging = 4);
	virtual ~SequentialGP();

	void computePosterior(const LikelihoodType& noiseModel);
	void computePosterior(const ivec& LikelihoodModel, const Vec<LikelihoodType *> noiseModels);
	void resetPosterior();
	void recomputePosterior();
	
	void computePosteriorFixedActiveSet(const LikelihoodType& noiseModel, ivec iActive);
	void recomputePosteriorFixedActiveSet(const LikelihoodType& noiseModel);
	
	/**
	 * Make predictions at a set of locations Xpred. The mean and variance
	 * are returned in the Mean and Variance vectors. To use a different  
	 * covariance function to the training one (useful to do non noisy predictions),
	 * you can pass an optional CovarianceFunction object.
	 */  
	void makePredictions(vec& Mean, vec& Variance, const mat& Xpred, CovarianceFunction &cf) const;
	void makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const;
	
	
	vec simulate(const mat& Xpred, bool approx) const;

	void setLikelihoodType(LikelihoodCalculation lc);

	vec getParametersVector() const;
	void setParametersVector(const vec p);

	double objective() const;
	vec gradient() const;
	void displayModelParameters() const;
	
	/**
	 * Accessors
	 */
	ivec getActiveSetIndices() { return idxActiveSet; }
	mat  getActiveSetLocations() { return ActiveSet; }
	
	void setActiveSetSize(int n) { maxActiveSet = n; }
	
	// This is used by ModelTrainer to update the model after
	// parameters have been changed.
	void updateModel();
	
	void setSelectiveSweep(bool b);
	
	
	
private:

	inline void addOne(int index, const LikelihoodType& noiseModel, const bool fixActiveSet);
	void addOne_siteRemoval(int index);
	void addOne_cavity(double sig0, double &mu, double &sigx, mat &KX, mat& Xmat);
	void addOne_updateGammaEhat(double &gamma, vec &eHat, const double sig0, const mat KX);
	void addOne_removeExtraPoints(const bool fixActiveSet);
	void addOne_removeCollapsedPoints();

	void updateSparse(mat& KX, vec& eHat, const double gamma, const double qtp1, const double rtp1, const double currentMean, const double currentVar, const double logEvidence, const int index);
	void updateFull(const mat& KX, vec& eHat, const double gamma, const double qtp1, const double rtp1, const double sig0, const double currentMean, const double currentVar, const double logEvidence, const int index);
	void updateEP();

	void stabiliseCoefficients(double& q, double& r, double cavityMean, double cavityVar, double upperTolerance, double lowerTolerance);
	vec scoreActivePoints(ScoringMethod sm);
	int findLeastInformativeActivePoint(ScoringMethod sm);
	void deleteActivePoint(int index);

	double compEvidence() const;
	double compEvidenceApproximate() const;
	double compEvidenceUpperBound() const;
	vec gradientEvidence() const;
	vec gradientEvidenceApproximate() const;
	vec gradientEvidenceUpperBound() const;

	mat computeCholesky(const mat& iM) const;
	mat computeInverseFromCholesky(const mat& C) const;

	mat KB;           // covariance between BV
	mat Q;            // inverse covariance between BV
	vec Alpha;        // alphas for calculating mean
	mat C;            // for calculating variance
	mat ActiveSet;    // Active set
	ivec idxActiveSet;// Active 

	mat projectionP;  // projection matrix
	vec lambdaP;      // lambda diagonal
	vec alphaP;       // alphaP
	vec logZ;         // logZ

	mat& Locations;
	vec& Observations;

	CovarianceFunction& covFunc;

	int sizeActiveSet;
	int maxActiveSet;
	// int sizeActiveSetOld;
	double epsilonTolerance;
	bool momentProjection;
	int iterChanging;
	int iterFixed;
	

	LikelihoodCalculation likelihoodType;
};

#endif /*SEQUENTIALGP_H_*/
