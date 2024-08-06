#include "spline.h"

spline::spline(Vector& nodesX, Vector& nodesY)
{
	this->nNodes = nodesX.getN();
	this->nodes = nodesX.copyData();

	this->polynomials = naturalSpline(nodesX, nodesY);
}

/*spline::spline(Vector& x, Vector& y, int flag)
{
	{
		switch (flag)
		{
		case 1:
			return naturalSpline(x, y);
		case 2:
			return periodicSpline(x, y);
		case 3:
			return notaKnotSpline(x, y);
		default:
			return notaKnotSpline(x, y);
		}
	}
}*/


polynomial* naturalSpline(Vector& nodesX, Vector& nodesY)
{
	int n, N, i, indexCon, indexFirDer, indexSecDer;
	float nodeX;

	n = nodesX.getN();
	N = 4 * (n - 1);

	Matrix A = Matrix(N);
	Matrix perm = Matrix(N);
	Vector c = Vector(N);

	//insert nodes
	i = 0;
	while (i < n) {
		nodeX = nodesX[i];
		A.setElement(i, 0, nodeX * nodeX * nodeX);
		A.setElement(i, 1, nodeX * nodeX);
		A.setElement(i, 2, nodeX);
		A.setElement(i, 3, 1.0);
		c.setElement(i, nodesY[i]);

		// insert continuity in internal nodes 
		indexCon = i + n;
		// s1(x1)
		A.setElement(indexCon, i,      nodeX * nodeX * nodeX);
		A.setElement(indexCon, i + 1,  nodeX * nodeX);
		A.setElement(indexCon, i + 2,  nodeX);
		A.setElement(indexCon, i + 3,  1);
		// -s2(x1)
		A.setElement(indexCon, i + 4, -nodeX * nodeX * nodeX);
		A.setElement(indexCon, i + 5, -nodeX * nodeX);
		A.setElement(indexCon, i + 6, -nodeX);
		A.setElement(indexCon, i + 7, -1);
		//
		c.setElement(indexCon, 0);

		// insert continuity of first derivative
		indexFirDer = indexCon + (n - 2);
		//s1'(x1)
		A.setElement(indexFirDer, i,      3 * nodeX * nodeX);
		A.setElement(indexFirDer, i + 1,  2 * nodeX);
		A.setElement(indexFirDer, i + 2,  1);
		//s2'(x1)
		A.setElement(indexFirDer, i,     -3 * nodeX * nodeX);
		A.setElement(indexFirDer, i + 1, -2 * nodeX);
		A.setElement(indexFirDer, i + 2, -1);
		//
		c.setElement(indexFirDer, 0);

		// insert continuity of second derivative
		indexSecDer = indexFirDer + (n - 2);
		//s1''(x1)
		A.setElement(indexSecDer, i,      6 * nodeX);
		A.setElement(indexSecDer, i + 1,  2);
		//s2''(x1)
		A.setElement(indexSecDer, i + 4, -6 * nodeX);
		A.setElement(indexSecDer, i + 5, -2);
		//
		c.setElement(indexSecDer, 0);
		i++;
	}
	
	// natural conditions 
	A.setElement(n - 1, 0, 6 * nodesX[0]);
	A.setElement(n - 1, 1, 2);
	c.setElement(n - 1, 0);
	A.setElement(n, n - 4, 6 * nodesX[n - 1]);
	A.setElement(n, n - 4, 2);
	c.setElement(n, 0);

	A.permutation(perm);
	c = perm * c;

	A.jacobi(c);
	
}