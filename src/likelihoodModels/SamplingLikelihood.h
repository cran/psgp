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

#ifndef SAMPLINGLIKELIHOOD_H_
#define SAMPLINGLIKELIHOOD_H_

#include <itpp/itbase.h>


#include "LikelihoodType.h"

using namespace itpp;

class SamplingLikelihood : public LikelihoodType
{
public:
	SamplingLikelihood();
	virtual ~SamplingLikelihood();

	void setSamplingParameters(int Samples, int Cycles);

	virtual double modelFunction(const double x) const = 0;
	vec modelFunction(const vec x) const;

	virtual double updateCoefficients(double& K1, double& K2, double Observation, double ModelMean, double ModelVariance) const = 0;


protected:
	int numberSamples;
	int numberCycles;

};
#endif /*SAMPLINGLIKELIHOOD_H_*/

