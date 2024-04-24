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

	Vector x = Vector(3);
	x.setElement(0, -2);
	x.setElement(1, 1);
	x.setElement(2, 3);
	//x.print();

	Vector y = Vector(3);
	y.setElement(0, 3);
	y.setElement(1, -7);
	y.setElement(2, -5);

	polynomial la = lagrangePoly(x, y);
	la.print();



	


}

