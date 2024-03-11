#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"


static float fx(float x) {
	return x * x + cos(2 * x) - 3;
}

int main() {
	
	/*
	* 
	* NEED TO IMPLEMENT A FUCKING EIGENVECTOR ALGORITHM, EIGENVALUES ARE OK
	* 
	  LAGRANGE POLYNOMIALS ARE UNFINISHED

	  INTEGRATION FORMULAS ARE BEHAVING VERY WEIRDLY WTF???????????

	  */

	Matrix A = Matrix(2);
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

	Vector x = A.gradient(b, x0, 1e-6, 1e3);
	Vector x2 = A.gradient(b);
	x.print();
	x2.print();



	
	


}

