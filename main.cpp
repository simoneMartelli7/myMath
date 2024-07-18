#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"
#include "polynomial.h"
#include "ODE1.h"
#include "ODE.h"
#include "functionalVector.h"

float fx(float u) {
	return 4.0 / (u * u + 1);
}

float fy(float u) {
	return cos(u);
}

int main() {
	
	/*
	* 
	* NEED TO IMPLEMENT A FUCKING EIGENVECTOR ALGORITHM, EIGENVALUES ARE OK

	  INTEGRATION FORMULAS ARE BEHAVING VERY WEIRDLY WTF???????????

	  */

	/*Matrix A = Matrix(2);
	Vector cauchy = Vector(2);
	float t0 = 0;

	cauchy.setElement(0, 0);
	cauchy.setElement(1, 1);
	A.setElement(0, 0);
	A.setElement(1, 1);
	A.setElement(2, -2);
	A.setElement(3, -3);

	std::vector <std::function<float(float)>> testF = { fx, fy };
	functionalVector testV = functionalVector(testF);

	ODE testForcing = ODE(A, cauchy, t0, testV);
	ODE testNoForcing = ODE(A, cauchy, t0);

	Vector* uEE = testForcing.explicitEuler(0.1, 10);
	Vector* uheun = testForcing.heun(0.1, 10);
	
	int i = 0;
	while (i < 5) {

		std::cout << "Euler: \n";
		uEE[i].print();
		std::cout << "Heun: \n";
		uheun[i].print();
		i++;
	}*/


	float integr = trapz(0, 1, 15, fx);

	std::cout << integr;
	
}

