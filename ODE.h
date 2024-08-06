#pragma once
#include "Matrix.h"
#include "functionalVector.h"

class ODE
{
private:
	int n;
	Matrix A;
	functionalVector g;
	Vector cauchy;
	float t0;
	int forcing; // set this value to 1 to evaluate the system with external forcing 
				 // specified by the functionalVector g

public: 
	ODE() {
		this->n = 2;
		this->A = Matrix(n);
		this->cauchy = Vector(n);
		this->t0 = 0;
		this->g = functionalVector();
		this->forcing = 0;

	}
	ODE(Matrix A, Vector cauchy, float t0) {
		this->n = A.getRows();
		this->A = A;
		this->cauchy = cauchy;
		this->t0 = t0;
		this->g = functionalVector();
		this->forcing = 0;
	}
	ODE(Matrix A, Vector cauchy, float t0, functionalVector g) {
		this->n = A.getRows();
		this->A = A;
		this->cauchy = cauchy;
		this->t0 = t0;
		this->g = g;
		this->forcing = 1;
	}

	// EULER'S METHOD 
	Vector* explicitEulerNoForcing(float deltaT, float T);
	Vector* explicitEulerForcing(float deltaT, float T);
	Vector* explicitEuler(float deltaT, float T) {
		return (forcing == 0 ? explicitEulerNoForcing(deltaT, T) : explicitEulerForcing(deltaT, T));
	}


	// HEUN'S METHOD           FIX THIS SHIT 
	Vector* heunNoForcing(float deltaT, float T);
	Vector* heunForcing(float deltaT, float T);
	Vector* heun(float deltaT, float T) {
		return (forcing == 0 ? heunNoForcing(deltaT, T) : heunForcing(deltaT, T));
	}

};

