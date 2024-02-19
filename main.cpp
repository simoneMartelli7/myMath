#include "Vector.h"
#include "Matrix.h"

int main() {
	
	/*float* data = new float[4];
	int j = 0;
	while (j < 4) {
		data[j] = -j;
		j++;
	}
	Vector u = Vector(4, data);
	u.print();
	Vector v = Vector(4);
	v.base(2);
	v.print();

	v = v * 2.1;

	float test = u.normInf();
	std::cout << u[3];*/

	Matrix A = Matrix(3);
	A.identity();
	A.setElement(2, 1, 0.3);
	Matrix B = Matrix(3);
	B.identity();
	B = B*3;
	B.setElement(2, 1.3);
	Matrix C = A * B;
	C.print();

	Vector u = Vector(3);
	u.base(1);
	u.setElement(2, 1);
	
	Matrix P = Matrix(3);
	std::vector<float> factors = C.gaussianElimination(P);
	C.print();


	
}