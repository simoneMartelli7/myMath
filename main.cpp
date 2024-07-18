#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"
#include "polynomial.h"
#include "ODE1.h"
#include "ODE.h"
#include "functionalVector.h"
#include "functionalMatrix.h"

float fx(float* u) {
	return 4.0 / (u[0] * u[1] + 1);
}

float fy(float* u) {
	return cos(u[0]) - sin(u[1]);
}

int main() {
	
	/*
	* 
	* NEED TO IMPLEMENT A FUCKING EIGENVECTOR ALGORITHM, EIGENVALUES ARE OK

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

	int i = 0;
	std::vector<std::function<float(float*)>> dummyF = { fx, fy, fy, fx };
	int n = 2;
	functionalMatrix test = functionalMatrix(n+1, n, dummyF);
	float* x = new float[2];
	x[0] = 0;
	x[1] = 3.141592;

	Matrix A = test.eval(x);
	A.print();
	

	
	
}

