#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"
#include "polynomial.h"
#include "ODE1.h"
#include "ODE.h"
#include "functionalVector.h"
#include "functionalMatrix.h"
#include "finiteDifferences.h"
#include "nonLinearSystems.h"
#include "spline.h"

float fx(float* u) {
	return 0;
}

float fy(float* u) {
	return log(u[0] + 1) + sin(u[0]);
}

float cosX(float x) {
	return cos(x);
}

int main() {
	
	/*
	* 
	* NEED TO IMPLEMENT A FUCKING EIGENVECTOR ALGORITHM, EIGENVALUES ARE OK
	* 
	* NON-LINEAR SYSTEMS NEED TO BE IMPROVED 
	* 
	* NEED TO ADDD MORE ODE METHODS, I'M NOT SURE HEUN IS WORKING CORRECTLY 

	  */
/*
	Matrix A = Matrix(2);
	Vector cauchy = Vector(2);
	float t0 = 0;

	cauchy.setElement(0, 1);
	cauchy.setElement(1, 2);
	A.setElement(0, 0);
	A.setElement(1, 1);
	A.setElement(2, -3);
	A.setElement(3, -3);

	std::vector <std::function<float(float*)>> testF = { fx, fy };
	functionalVector testV = functionalVector(testF);

	ODE testForcing = ODE(A, cauchy, t0, testV);
	ODE testNoForcing = ODE(A, cauchy, t0);

	Vector* uEE = testForcing.explicitEuler(0.1, 10);
	Vector* uheun = testForcing.heun(0.1, 10);
	
	int i = 0;

	while (i < 100) {

		uEE[i].save("euler.csv", 1);
		uheun[i].save("heun.csv", 1);
		i++;
	}*/

	/*int i = 0;
	std::vector<std::function<float(float*)>> dummyF = { fx, fy, fy, fx };
	int n = 2;
	functionalMatrix test = functionalMatrix(n+1, n, dummyF);
	float* x = new float[2];
	x[0] = 0;
	x[1] = 3.141592;

	Matrix A = test.eval(x);
	A.print();
	
	float df;
	Vector x0 = Vector(2);
	x0.setElement(0, 0);
	x0.setElement(1, 1);
	std::vector<std::function<float(float*)>> dummyF = { fx, fy };
	functionalVector testF = functionalVector(dummyF);

	//Vector test = newton(x0, testF, float(0.0001), 1000, float(0.001));
	//test.print();

	//Vector result = testF.eval(test.getData());
	//result.print();*/
	

	Vector nodesX = Vector(5);
	nodesX.setElement(0, 0.1);
	nodesX.setElement(1, 0.5);
	nodesX.setElement(2, 0.6);
	nodesX.setElement(3, 0.7);
	nodesX.setElement(4, 1.0);

	Vector nodesY = fillNodes(nodesX, cosX);


	nodesX.print();
	nodesY.print();
	
	spline test = spline(nodesX, nodesY, 0.0, centeredDifference(1.0, 0.001, cosX));

	std::cout << cosX(0.82) << " Interpolated value:  " << test.eval(0.82);

}

