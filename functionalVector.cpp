#include "functionalVector.h"

functionalVector::functionalVector(std::vector<std::function<float(float*)>> f)
{
	this->n = f.size();
	this->f = f;
}

Vector functionalVector::eval(float* x)
{
	Vector result = Vector(n);
	int i = 0;
	std::function<float(float*)> dummyF;

	while (i < n) {
		dummyF = getF(i);
		result.setElement(i, dummyF(x));
		i++;
	}

	return result;
}
