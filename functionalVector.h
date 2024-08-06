#pragma once
#include <vector>
#include <functional>
#include "Vector.h"
//for simplicity sake this accepts an array containing the variables of the functions 
//this is mainly because i can't be bothered to handle the case of a single variable, 
//care must be then given to the implementation of the functions themselves 
class functionalVector
{
private:
	int n;
	std::vector <std::function<float(float*)>> f;

public:
	functionalVector() {
		this->n = 2;
		std::vector<std::function<float(float*)>> dummy;
		this->f = dummy;
	}

	functionalVector(std::vector<std::function<float(float*)>> f);

	//returns a vector of the same dimension, the vector function is evaluated at point x
	Vector eval(float* x);

	std::function<float(float*)> getF(int i) { return f[i]; }
};

