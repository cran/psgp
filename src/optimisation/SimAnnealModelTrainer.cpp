#include "SimAnnealModelTrainer.h"

SimAnnealModelTrainer::SimAnnealModelTrainer(Optimisable& m) : ModelTrainer(m)
{
	algorithmName = "Simulated Annealing";
	stepSize = 0.02;
	coolingRate = 0.92;
	initialTemperature = 4.0;
	stepsPerTemperature = 1000;

}

SimAnnealModelTrainer::~SimAnnealModelTrainer()
{

}

void SimAnnealModelTrainer::Train(int numIterations){
	const double fmin = 1.0e-15;
	const double fmax = 1.0e100;
	vec dx, xnew, xsmallest, x = model.getParameters();
	int dims = x.size();
	vec xold = zeros(dims);
	double f = errorFunction(xold);
	double fnew = fmin, fold = fmax;
	double temperature, p, r, fsmallest = fmax;
		
	if(gradientCheck)
	{
		checkGradient();
	}
	
	for(int j = 1; j <= numIterations; j++)
	{
		temperature = initialTemperature * pow(coolingRate, (double)(j - 1.0));   

		for(int l = 1; l <= stepsPerTemperature; l++)
		{
			dx = stepSize * 2.0 * (randn(dims) - 0.5);

			xnew = x + dx;
			for(int i = 0; i < dims; i++)
			{
				if(xnew(i) <= eps)
				{
					xnew(i) = eps;
				}			
			}
 
			
			fnew = errorFunction(xnew);

			if(fnew < fsmallest)
			{
				fsmallest = fnew;
				xsmallest = xnew;
			}

			if((fnew >= fmin) && (fnew <= fmax))
			{
				if(fnew < f)
				{
					x = xnew;
					f = fnew;
				}
				else
				{
					p = exp((f - fnew) / temperature);
					r = randu();
					
					if(p >= r)
					{
						x = xnew;
						f = fnew;     
					}
				}
			}
		}
  
		if(display)
		{
			cout << "Cycle " << j;
			cout << "  Error " << fnew;
			cout << "  Temp. " << temperature;
			cout << endl;
		}


// do we need to check for some type of parameter tolerance values?
//		if((max(abs(x - xold)) < parameterTolerance) && (abs(f - fold) < errorTolerance))
//		{
//			functionValue = fsmallest;
//			model->setParameters(xsmallest);
//			return;
//		}
//		else
//		{
			xold = x;
			fold = f;
//		}
	} // for
		
//	if(display)
//	{
//		cout << "Warning: Maximum number of iterations has been exceeded" << endl;		
//	}

			functionValue = fsmallest;
			model.setParameters(xsmallest);

}


