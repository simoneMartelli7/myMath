#include "nonLinearSystems.h"
#include "finiteDifferences.h"

Vector newton(Vector x0, functionalVector F, float h, int maxIter, float tol)
{
	int i, j, n, k;
	float dfi_dxj, err;

	n = x0.getN();

	Matrix jacobian = Matrix(n);
	Vector distance = Vector(n);
	Vector f_x0 = Vector(n);

	err = 1;
	k = 0;
	while (k < maxIter && err > tol) {
		/*i = 0;
		while (i < n) {
			j = 0;
			std::function<float(float*)> fi = F.getF(i);
			while (j < n) {
				dfi_dxj = centeredDifferenceMulti(x0, j + 1, h, fi);
				jacobian.setElement(i, j, dfi_dxj);
				j++;
			}
			i++;
		}*/
		jacobian = Jacobian(F, h, x0);

		f_x0 = F.eval(x0.getData());
		f_x0 = f_x0 * (-1);
		distance = jacobian.qrSolve(f_x0); //HORRIBLE, BE BETTER
		err = distance.norm();
		x0 = x0 + distance;
		k++;
	}
	return x0;
}

Vector newtonMultiStep(int m, Vector x0, functionalVector F, float h, int maxIter, float tol)
{
	int n, k, i;
	float dfi_dxj, err;

	n = x0.getN();

	Matrix jacobian = Matrix(n);
	Vector yi = Vector(n);
	Vector yip1 = Vector(n);
	Vector f_x0 = Vector(n);

	k = 0;
	err = 1.0;

	while (k < maxIter && err > tol) {
		jacobian = Jacobian(F, h, x0); 

	}
	return Vector();
}


