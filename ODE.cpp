#include "ODE.h"

Vector* ODE::eulerExplicit(float deltaT, float T)
{
	float kTot = T / deltaT;
	Vector du;
	Vector* u = new Vector[kTot];
	int i;

	u[0] = cauchy;
	i = 1;
	while (i < kTot) {
		du = A * u[i - 1];
		du = du * deltaT;
		u[i] = u[i - 1] + du;
		i++;
	}

	return u;
}
