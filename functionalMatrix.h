#pragma once
#include <vector>
#include <functional>
#include "Matrix.h"

class functionalMatrix
{
private:
	int nRows;
	int nCols;
	// keep in mind that the program doesn't check for dimensional mismatch, be smart with it 
	std::vector<std::function<float(float*)>> F;
public:
	functionalMatrix() {
		this->nRows = 2;
		this->nCols = 2;
		std::vector<std::function<float(float*)>> dummy;
		this->F = dummy;
	}
	functionalMatrix(int nRows, int nCols, std::vector<std::function<float(float*)>> fGiven);

	// keep in mind that the program doesn't check for dimensional mismatch, be smart with it 
	Matrix eval(float* x);

	std::function<float(float*)> getF(int i) { return F[i]; }
};

