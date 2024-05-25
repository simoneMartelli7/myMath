#include "ODE1.h"

float* ODE1::explicitEuler(float T, float deltaT)
{
	float kTot = T / deltaT, du;
	float* u = new float[kTot];
	int i;

	u[0] = cauchy;

	i = 1;
	while (i < kTot) {
		du = f(t0 + deltaT * (i - 1), u[i - 1]);
		u[i] = u[i - 1] + deltaT * du;
		i++;
	}

	return u;
}
