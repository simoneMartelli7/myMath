#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"
#include "polynomial.h"
#include "ODE1.h"
#include "ODE.h"

float fx(float t, float u) {
	return sin(t) + u;
}

int main() {
	
	/*
	* 
	* NEED TO IMPLEMENT A FUCKING EIGENVECTOR ALGORITHM, EIGENVALUES ARE OK

	  INTEGRATION FORMULAS ARE BEHAVING VERY WEIRDLY WTF???????????

	  */

	Matrix A = Matrix(2);
	Vector cauchy = Vector(2);
	float t0 = 0;

	cauchy.setElement(0, 0);
	cauchy.setElement(1, 1);
	A.setElement(0, 0);
	A.setElement(1, 1);
	A.setElement(2, -2);
	A.setElement(3, -3);

	ODE test = ODE(A, cauchy, t0);

	Vector* uEE = test.eulerExplicit(0.1, 10);
	
	int i = 0;
	while (i < 25) {
		uEE[i].print();
		i++;
	}


}

