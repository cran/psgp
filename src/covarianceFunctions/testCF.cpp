#include <iostream>

#include "CovarianceFunction.h"
#include "GaussianCF.h"
#include <itpp/itbase.h>

using namespace std;
using namespace itpp;

int main()
{
	GaussianCF gcf(3.3, 4.5);

	
	cout << "starting" << endl;
	
	mat X("1 1;2 2;3 3;4 4;5 5");
	mat Z = zeros(5,5);

	gcf.computeSymmetric(Z, X);
	cout << Z << endl;

	gcf.getParameterPartialDerivative(Z, 0, X);
	cout << Z << endl;

	gcf.getParameterPartialDerivative(Z, 1, X);
	cout << Z << endl;
	
	cout << gcf.getParameters() << endl;
	gcf.displayCovarianceParameters();
	
}
