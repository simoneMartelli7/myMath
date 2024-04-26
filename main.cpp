#include "Vector.h"
#include "Matrix.h"
#include "approximation.h"
#include "polynomial.h"


static float fx(float x) {
	return sin(x);
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

	Vector x = createNodes(-0.1, 0.1, 6);
	x.print();

	Vector y = fillNodes(x, fx);

	polynomial la = lagrangePoly(x, y);
	la.print();




	


}

