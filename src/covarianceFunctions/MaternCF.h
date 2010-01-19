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

#ifndef MATERNCF_H_
#define MATERNCF_H_

#include "CovarianceFunction.h"

#include <cmath>
#include <cassert>
#include <itpp/itbase.h>


// the Matern covariance needs modified bessel functions
// if we're compiling in the R environment then we should use the built-in
// R functions, otherwise we will use the GSL.
#ifdef USING_R
#include <Rmath.h>
#else
#include <gsl/gsl_sf_bessel.h>
#endif

using namespace itpp;

class MaternCF : public CovarianceFunction
{
public:
	MaternCF(double lengthscale, double var, double kappa);
	MaternCF(vec parameters);
	
	virtual ~MaternCF();
	
	inline double computeElement(const vec& A, const vec& B) const;
	inline double computeDiagonalElement(const vec& A) const;
	void getParameterPartialDerivative(mat& PD, const int parameterNumber, const mat& X) const;
	
	void setParameter(int parameterNumber, const double value);
	double getParameter(int parameterNumber) const;
	
	string getParameterName(int parameterNumber) const;

	

	inline double computeBessel(const double x ,const double smoothness) const;
	mat besselkMat(const mat& X, double nu) const;
	
	inline double calcMatern(const vec& V) const;
	inline double calcGaussian(const vec& V) const;
	inline double calcMaternDiag() const;
	

	
	double variance;
	double range;
	double smoothness;
};

#endif /*MATERNCF_H_*/
