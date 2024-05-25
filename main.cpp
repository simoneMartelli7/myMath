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

	float cauchy = 0;
	float t0 = 15;

	ODE1 test = ODE1(fx, cauchy, t0);

	float T = 1;
	float deltaT = 0.1;

	float* u = test.explicitEuler(T, deltaT);
	float* err = new float[T / deltaT];

	int i = 0;
	while (i < T / deltaT) {
		float uAn = uEx(t0 + deltaT * i);
		err[i] = abs(uAn - u[i]) / abs(uAn);
		i++;
	}

	i = 0;
	while (i < T / deltaT) {
		std::cout << " " << u[i] << " ";
		i++;
	}
	


}

