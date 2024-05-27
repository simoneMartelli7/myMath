#pragma once

#include <functional>


class ODE1
{
private:
	std::function<float(float, float)> f;
	float cauchy;
	float t0;

public:
	ODE1() {
		this->f = 0;
		this->cauchy = 0;
		this->t0 = 0;
	}

	ODE1(std::function<float(float, float)> fGiven, float cauchy, float t0) {
		this->f = fGiven;
		this->cauchy = cauchy;
		this->t0 = t0;
	}

	// SETTER
	void setCauchy(float value) { this->cauchy = value; };
	void setT0(float value)     { this->t0 = value; };
	void changeF(std::function<float(float, float)> newF) { this->f = newF; };


	// EXPLIXCIT METHODS 
	float* explicitEuler(float T, float deltaT);
	float* heun(float T, float deltaT);
	float* modifiedEuler(float T, float deltaT);
	float* rungeKutta(float T, float deltaT, int s, float* a, float* b, float* c);

	// IMPLICIT METHODS 


};

