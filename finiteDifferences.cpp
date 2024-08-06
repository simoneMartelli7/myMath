#include "finiteDifferences.h"
#include "Vector.h"

float centeredDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0 + delta) - f(x0 - delta)) / delta;
}

float centeredDifference(std::function<float(float)> f, float xip, float xim)
{
	return (f(xip) - f(xim)) / (xip - xim);
}

float centeredDifference4(std::function<float(float)> f, float xip2, float xip, float xim, float xim2)
{
	return (-f(xip2) + 8 * f(xip) - 8 * f(xim) + f(xim2)) / (12 * (xip2 - xip));
}

float backwardDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0) - f(x0 - delta)) / delta;
}

float backwardDifference(std::function<float(float)> f, float xi, float xim)
{
	return (f(xi) - f(xim)) / (xi - xim);
}

float bacwardDifference3(std::function<float(float)> f, float xip, float xi, float xim, float xim2)
{
	return (2 * f(xip) + 3 * f(xi) - 6 * f(xim) + f(xim2)) / (6 * (xip - xi));
}

float forwardDifference(float x0, float delta, std::function<float(float)> f)
{
	return 0.5 * (f(x0 + delta) - f(x0)) / delta;
}

float forwardDifference(std::function<float(float)> f, float xi, float xip)
{
	return (f(xip) - f(xi)) / (xip - xi);
}

float forwardDifference3(std::function<float(float)> f, float xip2, float xip, float xi, float xim)
{
	return (-f(xip2) + 6 * f(xip) - 3 * f(xi) - 2 * f(xim)) / (6 * (xip - xi));
}

float centered2Difference(float x0, float delta, std::function<float(float)> f)
{
	return 0.25 * (f(x0 + delta) - 2 * f(x0) + f(x0 - delta)) / (delta * delta);
}

float centeredDifferenceMulti(Vector x0, int i, float h, std::function<float(float*)> f)
{
	Vector xp = Vector(x0.getN());
	Vector xm = Vector(x0.getN());
	xp.base(i);
	xp = xp * h;

	xm = x0 - xp;
	xp = x0 + xp;

	return 0.5 * (f(xp.getData()) - f(xm.getData())) / h;
}

Matrix Jacobian(functionalVector F, float h, Vector x0)
{
	int i, j, n;
	float dfi_dxj;

	n = x0.getN();
	Matrix jacobian = Matrix(n);
	i = 0;

	while (i < n) {
		j = 0;
		std::function<float(float*)> fi = F.getF(i);
		while (j < n) {
			dfi_dxj = centeredDifferenceMulti(x0, j + 1, h, fi);
			jacobian.setElement(i, j, dfi_dxj);
			j++;
		}
		i++;
	}

	return jacobian;
}




