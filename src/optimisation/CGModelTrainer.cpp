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
 // Conjugate Gradient algorithm based heavily on the Netlab
 // toolbox originally written in Matlab by Prof Ian Nabney
 //
 // History:
 //
 // 6 Dec 2008 - BRI - First implementation. 
 //

#include "CGModelTrainer.h"

CGModelTrainer::CGModelTrainer(Optimisable& m) : ModelTrainer(m)
{
	algorithmName = "Conjugate Gradient";
	lineMinimiserParameterTolerance = 1.0e-4;
}

CGModelTrainer::~CGModelTrainer()
{

}

void CGModelTrainer::Train(int numIterations){
	vec d, gradnew, gradold, xold, line_sd, x = getParameters();
	double br_min, br_max, gg, lmin, gamma, fold = 0.0, fnew = 0.0;
	
	if(gradientCheck)
	{
		checkGradient();
	}

	fnew = errorFunction(x);
	gradnew = errorGradients(x);
	d = -gradnew; // Initial search direction

	br_min = 0;
	br_max = 1.0;	// Initial value for maximum distance to search along

	for(int j = 1; j <= numIterations; j++)
	{
		xold = x;
		fold = fnew;
		gradold = gradnew;

		gg = dot(gradold, gradold);
		if (gg == 0.0)
		{
			functionValue = fnew;
			setParameters(x);
			return;
		}

		// This shouldn't occur, but rest of code depends on d being downhill
		if (dot(gradnew, d) > 0)
		{
			d = -d;
			if(display)
			{
				cout << "Warning: search direction uphill in conjugate gradient" << endl;		
			}
		}

		line_sd = d / sqrt(sum(elem_mult(d, d)));
		lineMinimiser(fnew, lmin, fold, xold, line_sd);
		x = xold + lmin * line_sd;

 		if((max(abs(x - xold)) < parameterTolerance) && (abs(fnew - fold) < errorTolerance))
		{
			functionValue = fnew;
			// RB
			setParameters(x);
			return;
		}

		gradnew = errorGradients(x);

		// Use Polak-Ribiere formula to update search direction
		// should we include an option to use other update formulas?
		gamma = dot((gradnew - gradold), gradnew) / gg;
		d = (d * gamma) - gradnew;

		if(display)
		{
			cout << "Cycle " << j;
			cout << "  Error " << functionValue << endl;
		}

	} // for
		
	if(display)
	{
		cout << "Warning: Maximum number of iterations has been exceeded" << endl;		
	}
	functionValue = fold;
	setParameters(x);
}


