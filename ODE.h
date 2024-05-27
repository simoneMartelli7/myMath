#pragma once
#include "Matrix.h"
class ODE
{
private:
	int n;
	Matrix A;
	Vector cauchy;
	float t0;
public: 
	ODE() {
		this->n = 2;
		this->A = Matrix(n);
		this->cauchy = Vector(n);
		this->t0 = 0;
	}
	ODE(Matrix A, Vector cauchy, float t0) {
		this->n = A.getRows();
		this->A = A;
		this->cauchy = cauchy;
		this->t0 = t0;
	}

	// EXPLICIT METHODS

	Vector* eulerExplicit(float deltaT, float T);
};

