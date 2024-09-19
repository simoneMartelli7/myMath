#include "spline.h"

// this is an helper function that compiles the system necessary to obtain the coefficients of the 
// various polynomials; works directly on the argument and fills only interpolation, continuity and 
// continuity of the derivatives rows, the two additional conditions must be filled indipendently for each case 
// also fills the c vector 
void fillMatrix(Matrix& A, Vector& c, Vector& x, Vector& y)
{
	int n, N, i, j;
	int indexFirstDer, indexSeconDer;
	float xValue;

	n = x.getN();
	N = c.getN();

	A.zero();
	c.zero();
	
	indexFirstDer = 2 * (n - 1) - 1;
	indexSeconDer = 3 * (n - 1) - 2;

	// first and last node interpolation is handled separately to facilitate the rest of the rows 
	xValue = x.getElement(0);
	A.setElement(0, 0, xValue * xValue * xValue);
	A.setElement(0, 1, xValue * xValue);
	A.setElement(0, 2, xValue);
	A.setElement(0, 3, 1);
	c.setElement(0, y.getElement(0));

	i = 2 * (n - 1) - 1;
	j = (n - 2) * (n - 1);
	xValue = x.getElement(n - 1);
	A.setElement(i, j, xValue * xValue * xValue);
	A.setElement(i, j + 1, xValue * xValue);
	A.setElement(i, j + 2, xValue);
	A.setElement(i, j + 3, 1);
	c.setElement(i, y.getElement(n - 1));

	i = 1;
	while (i < n - 1) {
		xValue = x.getElement(i);

		j = 2 * (i - 1) + 1;

		//interpolation and continuity
		A.setElement(j, 4 * (i - 1),     xValue * xValue * xValue);
		A.setElement(j + 1, 4 * i,       xValue * xValue * xValue);
		A.setElement(j, 4 * (i - 1) + 1, xValue * xValue);
		A.setElement(j + 1, 4 * i + 1,   xValue * xValue);
		A.setElement(j, 4 * (i - 1) + 2, xValue);
		A.setElement(j + 1, 4 * i + 2,   xValue);
		A.setElement(j, 4 * (i - 1) + 3, 1);
		A.setElement(j + 1, 4 * i + 3,   1);

		c.setElement(j,     y.getElement(i));
		c.setElement(j + 1, y.getElement(i));


		// continuity of first derivative
		A.setElement(indexFirstDer + i, 4 * (i - 1),         3 * xValue * xValue);
		A.setElement(indexFirstDer + i, 4 * i,           -3 * xValue * xValue);
		A.setElement(indexFirstDer + i, 4 * (i - 1) + 1,     2 * xValue);
		A.setElement(indexFirstDer + i, 4 * i + 1,      -2 * xValue);
		A.setElement(indexFirstDer + i, 4 * (i - 1) + 2,     1);
		A.setElement(indexFirstDer + i, 4 * i + 2,      -1);
		
		c.setElement(indexFirstDer + i, 0);
		c.setElement(indexFirstDer + i + 1, 0);


		//continuity of second derivative
		A.setElement(indexSeconDer + i, 4 * (i - 1),     6 * xValue);
		A.setElement(indexSeconDer + i, 4 * i,      -6 * xValue);
		A.setElement(indexSeconDer + i, 4 * (i - 1) + 1,     2);
		A.setElement(indexSeconDer + i, 4 * i + 1,      -2);

		c.setElement(indexSeconDer + i, 0);


		i++;
	}
}

polynomial* naturalSpline(Vector& nodesX, Vector& nodesY)
{
	int n, N, i, indexCon, indexFirDer, indexSecDer;
	float nodeX;

	n = nodesX.getN();
	N = 4 * (n - 1);

	Matrix A = Matrix(N);
	Matrix perm = Matrix(N);
	Vector c = Vector(N);
	Vector coeffs;
	polynomial localPoly = polynomial(3);
	float* dummy = new float [4] {0.0};
	float* interval = new float[2] {0.0};
	polynomial* result = new polynomial[n];

	
	fillMatrix(A, c, nodesX, nodesY);

	// natural conditions 
	A.setElement(N - 2, 0, 6 * nodesX[0]);
	A.setElement(N - 2, 1, 2);
	c.setElement(N - 2, 0);
	A.setElement(N - 1, N - 4, 6 * nodesX[n - 1]);
	A.setElement(N - 1, N - 3, 2);
	c.setElement(N - 1, 0);


	coeffs = A.paluSolve(c);

	i = 0;
	while (i < n) {
		dummy[0] = coeffs[4 * i + 3];
		dummy[1] = coeffs[4 * i + 2];
		dummy[2] = coeffs[4 * i + 1];
		dummy[3] = coeffs[4 * i];
		localPoly = polynomial(3, dummy);
		result[i] = localPoly;
		i++;
	}
	return result;
}

polynomial* periodicSpline(Vector& nodesX, Vector& nodesY)
{
	int n, N, i, indexCon, indexFirDer, indexSecDer;
	float nodeX;

	n = nodesX.getN();
	N = 4 * (n - 1);

	Matrix A = Matrix(N);
	Matrix perm = Matrix(N);
	Vector c = Vector(N);
	Vector coeffs;
	polynomial localPoly = polynomial(3);
	polynomial* result = new polynomial[n];
	float* dummy = new float [4] {0.0};

	fillMatrix(A, c, nodesX, nodesY);

	//periodic conditions 
	// s'(x1) = s'(xn)
	//s''(x1) = s''(xn)
	nodeX = nodesX[0];
	A.setElement(N - 2, 0, 3 * nodeX * nodeX);
	A.setElement(N - 2, 1, 2 * nodeX);
	A.setElement(N - 2, 2, 1);
	A.setElement(N - 1, 0, 6 * nodeX);
	A.setElement(N - 1, 1, 2);

	nodeX = nodesX[n];
	A.setElement(N - 2, N - 4, -3 * nodeX * nodeX);
	A.setElement(N - 2, N - 3, -2 * nodeX);
	A.setElement(N - 2, N - 2, -1);
	A.setElement(N - 1, N - 4, -6 * nodeX);
	A.setElement(N - 1, N - 3, -2);
	
	
	c.setElement(N - 1, 0);
	c.setElement(N - 2, 0);

	coeffs = A.paluSolve(c);

	i = 0;
	while (i < n) {
		dummy[0] = coeffs[4 * i + 3];
		dummy[1] = coeffs[4 * i + 2];
		dummy[2] = coeffs[4 * i + 1];
		dummy[3] = coeffs[4 * i];
		localPoly = polynomial(3, dummy);
		result[i] = localPoly;
		i++;
	}
	return result;
}

polynomial* notaKnotSpline(Vector& nodesX, Vector& nodesY)
{
	int n, N, i, indexCon, indexFirDer, indexSecDer;
	float nodeX;

	n = nodesX.getN();
	N = 4 * (n - 1);

	Matrix A = Matrix(N);
	Matrix perm = Matrix(N);
	Vector c = Vector(N);
	Vector coeffs;
	polynomial localPoly = polynomial(3);
	polynomial* result = new polynomial[n - 1];
	float* dummy = new float [4] {0.0};

	fillMatrix(A, c, nodesX, nodesY);


	// not a knot conditions
	// s_0'''(x1) = s_1'''(x1) etc.
	A.setElement(N - 2, 0, 6);
	A.setElement(N - 2, 4, -6);
	c.setElement(N - 2, 0);
	A.setElement(N - 1, N - 8, 6);
	A.setElement(N - 1, N - 4, -6);
	c.setElement(N - 1, 0);

	coeffs = A.paluSolve(c);

	coeffs.print();

	i = 0;
	while (i < n - 1) {
		dummy[0] = coeffs[4 * i + 3];
		dummy[1] = coeffs[4 * i + 2];
		dummy[2] = coeffs[4 * i + 1];
		dummy[3] = coeffs[4 * i];
		localPoly = polynomial(3, dummy);
		result[i] = localPoly;
		i++;
	}
	return result;
}

polynomial* clampedSpline(Vector& nodesX, Vector& nodesY, float y10, float y1n)
{

	int n, N, i, indexCon, indexFirDer, indexSecDer;
	float nodeX;

	n = nodesX.getN();
	N = 4 * (n - 1);

	Matrix A = Matrix(N);
	Matrix perm = Matrix(N);
	Vector c = Vector(N);
	Vector coeffs;
	polynomial localPoly = polynomial(3);
	polynomial* result = new polynomial[n];
	float* dummy = new float [4] {0.0};

	fillMatrix(A, c, nodesX, nodesY);

	//clamped condition
	// s'(x1) = y0n
	//s'(xn) = y1n
	nodeX = nodesX[0];
	A.setElement(N - 1, 0, 3 * nodeX * nodeX);
	A.setElement(N - 1, 1, 2 * nodeX);
	A.setElement(N - 1, 2, 1);
	c.setElement(N - 1, y10);

	nodeX = nodesX[n];
	A.setElement(N - 2, N - 4, 3 * nodeX * nodeX);
	A.setElement(N - 2, N - 3, 2 * nodeX);
	A.setElement(N - 2, N - 2, 1);
	c.setElement(N - 2, y1n);


	coeffs = A.paluSolve(c);

	i = 0;
	while (i < n) {
		dummy[0] = coeffs[4 * i + 3];
		dummy[1] = coeffs[4 * i + 2];
		dummy[2] = coeffs[4 * i + 1];
		dummy[3] = coeffs[4 * i];
		localPoly = polynomial(3, dummy);
		result[i] = localPoly;
		i++;
	}
	return result;
}




spline::spline(Vector& nodesX, Vector& nodesY)
{
	this->nodes = nodesX;
	this->polynomials = notaKnotSpline(nodesX, nodesY);
}

spline::spline(Vector& nodesX, Vector& nodesY, int flag)
{
	{
		switch (flag)
		{
		case 1:
			this->nodes = nodesX;
			this->polynomials = naturalSpline(nodesX, nodesY);
			break;
		case 2:
			this->nodes = nodesX;
			this->polynomials = periodicSpline(nodesX, nodesY);
			break;
		case 3:
			this->nodes = nodesX;
			this->polynomials = notaKnotSpline(nodesX, nodesY);
			break;
		default:
			this->nodes = nodesX;
			this->polynomials = naturalSpline(nodesX, nodesY);
			break;
		}
	}
}

spline::spline(Vector& nodesX, Vector& nodesY, float y10, float y1n)
{
	this->nodes = nodesX;
	this->polynomials = clampedSpline(nodesX, nodesY, y10, y1n);
}

int spline::findInterval(float x0)
{
	float a, b;
	int index, flag, n;

	n = nodes.getN();
	a = nodes[0];
	b = nodes[n - 1];

	index = int(n / 2);
	flag = 1;

	if (x0 < a) {
		std::cout << "The extremes of interpolation are [" << a << ", " << b << "]";
		return 0;
	}

	if (x0 > b) {
		std::cout << "The extremes of interpolation are [" << a << ", " << b << "]";
		return n - 2;
	}

	while (flag) {
		a = nodes[index];
		b = nodes[index + 1];

		if (a < x0 && x0 < b) {
			flag = 0;
			return index;
		}
		else if (x0 < a) {
			index = int(index / 2);
		}
		else {
			index = index + int(index / 2);
		}
	}
}

float spline::eval(float x0)
{
	int index = findInterval(x0);

	return polynomials[index].eval(x0);
}




