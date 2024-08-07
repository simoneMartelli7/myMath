#include "approximation.h"
#include "finiteDifferences.h"

polynomial polynomialFitting(Vector& x, Vector& y)
{
	int i, j;
	int deg = x.getN();
	Matrix A = Matrix(deg);
	Vector C = Vector(deg);

	if (deg != y.getN()) {
		std::cerr << "Error: mismatch in the nodes";
		exit(-1);
	}

	i = 0;
	while (i < deg) {
		A.setElement(i, 0, 1);
		j = 1;
		while (j < deg) {
			A.setElement(i, j, pow(x.getElement(i), j));
			j++;
		}
		i++;
	}

	C = A.paluSolve(y);

	return polynomial(deg - 1, C.copyData());
}

float polynomialFitting(Vector& x, Vector& y, float x0)
{
	polynomial poly = polynomialFitting(x, y);

	return poly.eval(x0);
}

polynomial leastSquareFitting(Vector& x, Vector& y, int deg)
{
	int i, j;
	Matrix A = Matrix(x.getN(), deg+1);
	Vector C = Vector(deg+1);

	i = 0;
	while (i < x.getN()) {
		A.setElement(i, 0, 1);
		j = 1;
		while (j < deg+1) {
			A.setElement(i, j, pow(x.getElement(i), j));
			j++;
		}
		i++;
	}

	Matrix At = A.transpose();
	A = At * A;
	y = At * y;

	C = A.paluSolve(y);

	return polynomial(deg, C.copyData());
}


polynomial lagrangePoly(Vector& x, Vector& y)
{
	float denominator, x_k;
	int i, k;
	polynomial result = polynomial();
	polynomial dummy1 = polynomial(1);

	result.clear();
	k = 0;
	while (k < x.getN()) {
		x_k = x[k];
		i = 0;
		polynomial dummy2 = polynomial(0);
		dummy2.setCoeff(0, 1);
		denominator = 1;
		while (i < k) {
			dummy1.setCoeff(0, -x[i]);
			dummy1.setCoeff(1, 1);
			denominator = denominator * (x_k - x[i]);
			i++;
			dummy2 = dummy2.product(dummy1);
		}
		i++;
		while (i < x.getN()) {
			dummy1.setCoeff(0, -x[i]);
			dummy1.setCoeff(1, 1);
			denominator = denominator * (x_k - x[i]);
			i++;
			dummy2 = dummy2.product(dummy1);
		}
		dummy2 = dummy2.scalarProduct(y[k] / denominator);
		result = result.sum(dummy2);
		k++;
	}

	return result;
	
}


//these are mainly for debugging
Vector createNodes(float lowerLimit, float upperLimit, int nIntervals)
{
	float delta = (upperLimit - lowerLimit) / (nIntervals-1);
	int i = 0;
	Vector nodes = Vector(nIntervals);

	while (i < nIntervals) {
		nodes.setElement(i, lowerLimit + delta * i);
		i++;
	}

	return nodes;
}

Vector fillNodes(Vector& x, std::function<float(float)> f)
{
	int i, n = x.getN();
	Vector y = Vector(n);

	i = 0;
	while (i < n) {
		y.setElement(i, f(x[i]));
		i++;
	}

	return y;
}


float trapz(Vector& x, Vector& y)
{
	int i, n = x.getN();
	float integral = 0;

	if (n != y.getN()) {
		std::cerr << "Input Error: incompatible dimensions";
		exit(-1);
	}

	i = 1;
	while (i < n) {
		integral += 0.5 * ((x[i] - x[i - 1]) * (y[i] + y[i - 1]));
		i++;
	}

	return integral;
}

float trapz(float lowerLimit, float upperLimit, int N, std::function<float(float)> f)
{
	float h, integral;
	int i;

	h = (upperLimit - lowerLimit) / N;

	i = 1;
	integral = f(upperLimit) + f(lowerLimit);

	while (i < N) {
		integral += 2*f(lowerLimit + h * i);
		i++;
	}

	return 0.5 * integral * h;
}

float simpson(Vector& x, Vector& y)
{
	float integral, delta, k1, k2;
	int i, n = x.getN();

	if (n != y.getN()) {
		std::cerr << "Input Error: incompatible dimensions";
		exit(-1);
	}

	integral = 0;
	delta = x[1] - x[0];
	k1 = 1.0 / 6.0 * delta;
	k2 = 4.0 * k1;

	i = 2;
	while (i < x.getN()) {
		integral += (k1 * (y[i - 2] + y[i]) + k2 * y[i - 1]);
		i++;
	}

	return integral;
}

float simpson(float lowerLimit, float upperLimit, int N, std::function<float(float)> f)
{
	float integral, h;
	int k;

	integral = f(lowerLimit) + f(upperLimit);
	h = (upperLimit - lowerLimit) / N;

	integral += 4*f(lowerLimit + 0.5 * h);
	k = 1;
	while (k < N) {
		integral += (2 * f(lowerLimit + k * h) + 4 * f(lowerLimit + (k + 0.5 ) * h));
		k++;
	}

	return integral * h / 6;

}

float bisection(float a, float b, std::function<float(float)> fx, float tol)
{
	int k = ceil(log2((b - a) / tol));
	int i;
	float xi, fa, fb, fxi;

	fa = fx(a);
	fb = fx(b);

	if (fa * fb > 0) {
		std::cerr << "Invalid interval";
		exit(-1);
	}

	xi = 0.5 * (a + b);
	fxi = fx(xi);

	i = 0;
	while (i < k && abs(fxi) > tol) {
		xi = 0.5 * (a + b);
		fxi = fx(xi);
		(fa * fxi < 0) ? b = xi : a = xi;
		
		i++;
	}

	return xi;
}

float newtonDiscrete(float x0, std::function<float(float)> f, float tol, int maxIter)
{
	float fx, df, delta, x;
	int k;

	fx = f(x0);
	delta = 1e-3;
	df = centeredDifference(x0, delta, f);
	x = x0 - fx / df;
	
	k = 0;
	while (k < maxIter && abs(fx) > tol) {
		x0 = x;
		fx = f(x0);
		df = centeredDifference(x0, delta, f);
		x = x0 - fx / df;
		k++;
	}
	std::cout << "N. of iterations: " << k;

	return x;
}

float secant(float x1, float x2, std::function<float(float)> f, float tol, int maxIter)
{
	float fx1, fx2, x;
	int k;

	fx1 = f(x1);
	fx2 = f(x2);
	
	k = 0;
	while (k < maxIter && abs(fx2) > tol) {
		x = x2 - fx2 * (x2 - x1) / (fx2 - fx1);
		x1 = x2;
		x2 = x;
		fx1 = fx2;
		fx2 = f(x2);
		k++;
	}

	return x;
}
