/***************************************************************************
 *   AstonGeostats, algorithms for low-rank geostatistical models          *
 *                                                                         *
 *   Copyright (C) Ben Ingram, 2008                                        *
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
#include "optimisation/Optimisable.h"
#include "covarianceFunctions/CovarianceFunction.h"
#include "likelihoodModels/LikelihoodType.h"
#include "likelihoodModels/GaussianSampLikelihood.h"

#include <cassert>

class SequentialGP : public ForwardModel, public Optimisable
{
public:
	SequentialGP(int Inputs, int Outputs, mat& Xdata, vec& ydata, CovarianceFunction& cf);
	virtual ~SequentialGP();

	void computePosterior(const LikelihoodType& noiseModel);
	void computePosterior(const ivec& LikelihoodModel, const Vec<LikelihoodType *> noiseModels);
	void makePredictions(vec& Mean, vec& Variance, const mat& Xpred) const;

	vec getParametersVector() const;
	void setParametersVector(const vec p);

	double objective() const;
	vec gradient() const;

	void estimateParameters();

	void displayCovarianceParameters() const;

private:

	void addOne(const vec& X, const double Observation, const LikelihoodType& noiseModel);

	void sparseUpdate(mat& KX, mat& eHat, const double qtp1, const double rtp1);
	void fullUpdate(const mat& X, const mat& KX, mat& eHat, const double gamma, const double qtp1, const double rtp1, const double sig0);

	mat computeCholesky(const mat& iM) const;
	mat computeInverseFromCholesky(const mat& C) const;


	mat KB; // covariance between BV
	mat Q; // inverse covariance between BV
	mat Alpha; // alphas for calculating mean
	mat C; // for calculating variance
	mat ActiveSet; // Active set

	CovarianceFunction& covFunc;
	mat& Locations;
	vec& Observations;

	int sizeActiveSet;
	int maxActiveSet;
	double epsilonTolerance;

};

#endif /*SEQUENTIALGP_H_*/
