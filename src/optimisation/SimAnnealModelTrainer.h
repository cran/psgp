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

#ifndef SIMANNEALMODELTRAINER_H_
#define SIMANNEALMODELTRAINER_H_

#include <iostream>
#include <vector>
#include <string>
#include <itpp/itbase.h>

#include "Optimisable.h"
#include "ModelTrainer.h"


using namespace std;
using namespace itpp;

class SimAnnealModelTrainer : public ModelTrainer
{
public:
	SimAnnealModelTrainer(Optimisable& m);
	virtual ~SimAnnealModelTrainer();

	void Train(int numIterations);

	void setCoolingRate(double d) {coolingRate = d;};
	void setInitialTemperature(double d) {initialTemperature = d;};
	void setStepsPerTemperature(int i) {stepsPerTemperature = i;};
	void setStepSize(double d) {stepSize = d;};

protected:

private:
	double stepSize;
	double coolingRate;
	double initialTemperature;
	int stepsPerTemperature;

};


#endif /*SIMANNEALMODELTRAINER_H_*/
