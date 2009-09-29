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
 //
 // Gradient Descent algorithm based heavily on the Netlab
 // toolbox originally written in Matlab by Prof Ian Nabney
 //
 // History:
 //
 // 6 Dec 2008 - BRI - First implementation. 
 //

#include "GDModelTrainer.h"

GDModelTrainer::GDModelTrainer(Optimisable& m) : ModelTrainer(m)
{
	algorithmName = "Gradient Descent";
	eta = 0.01;
	mu = 0.1;
}

GDModelTrainer::~GDModelTrainer()
{

}

void GDModelTrainer::Train(int numIterations){
	// hardwired first
	bool computeObjectiveEachIteration = true;
	bool useLineMinimisation = true;
	
	if(useLineMinimisation)
	{
		algorithmName = "Gradient Descent (using line minimiser)";
	}
	
	vec x = getParameters();
	vec xold = x;
	vec dxold = zeros(x.size());
	vec grad;
	vec dx;
	double fold = 0.0;
	double fnew = 0.0;
	
	if(gradientCheck)
	{
		checkGradient();
	}
	
	if(computeObjectiveEachIteration)
	{
		fnew = errorFunction(x);
		fold = fnew;
	}

	for(int j = 1; j <= numIterations; j++)
	{
		xold = x;
		grad = errorGradients(x);
	
		if(useLineMinimisation)
		{
			double lmin;
			vec sd = - grad / sqrt(sum(elem_mult(grad, grad))); //% New search direction.
			fold = fnew;
			lineMinimiser(fnew, lmin, fold, x, sd);
			x = xold + (lmin * sd);	
		}
		else
		{
			if(true) // adaptable learning rate....
			{
	      	// Let learning rate decrease as 1/t
				mu = eta / j;
				mu = 0.5;
			}    		
			dx = (mu * dxold) - (eta * grad);
			x = x + dx;
			dxold = dx;
			if(computeObjectiveEachIteration)
			{
				fold = fnew;
				fnew = errorFunction(x);
			} // if
		} // if
					
		if(display)
		{
			cout << "Cycle " << j;
			cout << "  Error " << fnew;
			
			if(!useLineMinimisation)
			{
				cout << "  mu " << mu;
			}
			cout << endl;
		}

		if((max(abs(x - xold)) < parameterTolerance) && (abs(fnew - fold) < errorTolerance))
		{
			functionValue = fnew;
			return;
		}

	} // for

	if(computeObjectiveEachIteration)
	{
		functionValue = fnew;
	}
	else
	{
		fnew = errorFunction(x);
		functionValue = fnew;
	}
		
	if(display)
	{
		cout << "Warning: Maximum number of iterations has been exceeded" << endl;		
	}

}


