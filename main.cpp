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
	int i = 0;
	A.setElement(0, 3);
	A.setElement(1, 0);
	A.setElement(2, 1);
	A.setElement(3, 0);
	A.setElement(4, 3);
	A.setElement(5, 0);
	A.setElement(6, 1);
	A.setElement(7, 0);
	A.setElement(8, 3);


	A.print();
	Vector u = Vector(3);
	u.setElement(0, 1);
	u.setElement(1, 3);
	u.setElement(2, 3);
	u.print();
/*
	A = Matrix(2);
	A.setElement(0, 5);
	A.setElement(1, 4);
	A.setElement(2, 1);
	A.setElement(3, 2);

	float* eigenValues = new float[2] {0.0};

	std::unordered_map<float, Vector> eigen = A.qrEigen(10000, 1e-7, eigenValues);

	std::cout << eigenValues[0] << "\n";
	eigen[eigenValues[0]].print();

	std::cout << eigenValues[1] << "\n";
	eigen[eigenValues[1]].print();*/

	Vector x0 = Vector(3);
	Vector x = A.gradient(u, x0, 1e-7, 1e4);
	x.print();

}