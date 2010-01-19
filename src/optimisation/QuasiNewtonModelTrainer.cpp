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
 // Quasi-Newton algorithm based heavily on the Netlab
 // toolbox originally written in Matlab by Prof Ian Nabney
 //
 // History:
 //
 // 6 Dec 2008 - BRI - First implementation. 
 //
 
#include "QuasiNewtonModelTrainer.h"

QuasiNewtonModelTrainer::QuasiNewtonModelTrainer(Optimisable& m) : ModelTrainer(m)
{
	algorithmName = "Quasi-Newton";
	lineMinimiserParameterTolerance = 1.0e-2;	
}

QuasiNewtonModelTrainer::~QuasiNewtonModelTrainer()
{

}
void QuasiNewtonModelTrainer::Train(int numIterations){
	
	vec x = getParameters();
	vec d, gradold, xold, v, p, Gv, u;
	double lmin, vdotp, vGv, fold = 0.0;

	double fnew = errorFunction(x);
	vec gradnew = errorGradients(x);



	mat hessinv = eye(gradnew.size());
	double min_frac_change = 1.0e-4;
	
	if(gradientCheck)
	{
		checkGradient();
	}

	p = -gradnew; // Initial search direction

	for(int j = 1; j <= numIterations; j++)
	{
		xold = x;
		fold = fnew;
		gradold = gradnew;
		x = xold + p;
		fnew = errorFunction(x);

		if (dot(gradnew, p) >= 0)
		{
			p = -p;
			if(display)
			{
				cout << "Warning: search direction uphill in conjugate gradient" << endl;		
			}
		}

		if (fnew >= (fold + min_frac_change * dot(gradnew, p)))
		{
			lineMinimiser(fnew, lmin, fold, xold, p);
			x = xold + lmin * p;
			p = x - xold;
			fnew = functionValue;
		}


 		if((max(abs(x - xold)) < parameterTolerance) & (abs(fnew - fold) < errorTolerance))
		{
			functionValue = fnew;
			setParameters(x);
			return;
		}
		gradnew = errorGradients(x);

		v = gradnew - gradold;
		vdotp = dot(v, p);

		if ((vdotp*vdotp) > (eps*sum(elem_mult(v, v))*sum(elem_mult(p, p)))) 
		{
			Gv = hessinv*v;
			vGv = dot(v, Gv);
			u = (p / vdotp) - (Gv / vGv);
			// Use BFGS update rule
			// she we consider some other options for update rules?
			hessinv = hessinv + (outer_product(p, p) / vdotp) - (outer_product(Gv, Gv) / vGv) + (vGv * outer_product(u, u));
		}
		p = -(hessinv * gradnew);

		if(display)
		{
			cout << "Cycle " << j;
			cout << "  Error " << fnew << endl;
		}
	} // for
		
	if(display)
	{
		cout << "Warning: Maximum number of iterations has been exceeded" << endl;		
	}
	functionValue = fold;
	setParameters(x);
}


