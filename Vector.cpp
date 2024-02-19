#include "Vector.h"

Vector Vector::addition(Vector& v)
{
	if (v.getN() != getN()) {
		std::cerr << "Incompatible dimensions\n";
		exit(-1);
	}

	Vector result = Vector(v.getN());
	int j = 0;
	while (j < v.getN()) {
		result.setElement(j, getElement(j) + v.getElement(j));
		j++;
	}
	return result;
}

Vector Vector::scalarProduct(float alpha)
{
	Vector result = Vector(getN());
	int j = 0;
	while (j < n) {
		result.setElement(j, alpha * getElement(j));
		j++;
	}
	return result;
}

float Vector::dotProduct(Vector& v)
{
	float result = 0;
	int j = 0;
	while (j < n) {
		result += v.getElement(j) * getElement(j);
		j++;
	}
	return result;
}

float Vector::norm()
{
	int j = 0;
	float result = 0;
	while (j < n) {
		result += getElement(j) * getElement(j);
		j++;
	}
	return sqrt(result);
}

float Vector::norm1()
{
	int j = 0;
	float result = 0;
	while (j < n) {
		result += getElement(j);
		j++;
	}
	return result / n;
}

float Vector::normInf()
{
	float max = getElement(0);
	int j = 0;
	while (j < n) {
		j++;
		if (max - getElement(j) < 0) {
			max = getElement(j);
		}
	}
	return max;
}
