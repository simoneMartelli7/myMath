#pragma once
#include <vector>
#include <functional>
#include "Vector.h"

class functionalVector
{
private:
	int n;
	std::vector <std::function<float(float)>> f;

public:
	functionalVector() {
		this->n = 2;
		std::vector<std::function<float(float)>> dummy;
		this->f = dummy;
	}

	functionalVector(std::vector<std::function<float(float)>> f);

	//returns a vector whose of the same dimension, the vector function is evaluated at point x
	Vector eval(float x);

	std::function<float(float)> getF(int i) { return f[i]; }
};

