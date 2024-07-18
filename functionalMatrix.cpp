#include "functionalMatrix.h"

functionalMatrix::functionalMatrix(int nRows, int nCols, std::vector<std::function<float(float*)>> fGiven)
{
	this->nRows = nRows;
	this->nCols = nCols;

	this->F = fGiven;
}

Matrix functionalMatrix::eval(float* x)
{
	Matrix result = Matrix(nRows, nCols);
	int i = 0;
	std::function<float(float*)> dummyF;

	while (i < nRows * nCols) {
		dummyF = getF(i);
		result.setElement(i, dummyF(x));
		i++;
	}

	return result;
}
