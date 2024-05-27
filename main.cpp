#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"
#include "polynomial.h"
#include "ODE1.h"

float fx(float t, float u) {
	return sin(t) + u;
}

float uEx(float t) {
	return 0.5 * (exp(t) - sin(t) - cos(t));
}

int main() {
	
	/*
	* 
	* NEED TO IMPLEMENT A FUCKING EIGENVECTOR ALGORITHM, EIGENVALUES ARE OK

	  INTEGRATION FORMULAS ARE BEHAVING VERY WEIRDLY WTF???????????

	  */

	/*Matrix A = Matrix(2);
	int i = 0;
	A.setElement(0, 4);
	A.setElement(1, 1);
	A.setElement(2, 1);
	A.setElement(3, 3);

	Vector b = Vector(2);
	b.setElement(0, 1);
	b.setElement(1, 2);

	Vector x0 = Vector(2);
	x0.setElement(0, 2);
	x0.setElement(1, 1);

	Vector x = createNodes(-0.1, 0.1, 6);
	x.print();

	Vector y = fillNodes(x, fx);

	polynomial la = polynomial(3);
	la.setCoeff(3, 1);
	la.setCoeff(1, -1);
	la.setCoeff(0, -24);
	la.print();

	polynomial r = polynomial(1);
	r.setCoeff(0, -3);

	polynomial div = la.division(r);
	div.print();*/

	float cauchy = 1;
	float t0 = 0;

	ODE1 test = ODE1(fx, cauchy, t0);

	float T = 1;
	float deltaT = 0.1;

	float* uEE = test.explicitEuler(T, deltaT);
	
	int s = 3;
	float* a = new float[3];
	a[0] = 0.25;
	a[1] = 0;
	a[2] = 0.75;
	float* b = new float[3];
	b[0] = 0;
	b[1] = 1.0/3;
	b[2] = 2.0/3;
	float* c = new float[9] {0.0};
	c[3] = 1.0 / 3;
	c[6] = 0;
	c[7] = 2.0 / 3;

	float* uRK = test.rungeKutta(T, deltaT, 3, a, b, c);
	

	int i = 0;
	while (i < T / deltaT) {
		std::cout << " " << uEE[i] << " ";
		i++;
	}

	std::cout << "\n";
	i = 0;
	while (i < T / deltaT) {
		std::cout << " " << uRK[i] << " ";
		i++;
	}

	
	


}

