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

float* ODE1::heun(float T, float deltaT)
{
	float kTot = T / deltaT, du;
	float* u = new float[kTot];
	int i;

	u[0] = cauchy;

	i = 1;
	while (i < kTot) {
		du = f(t0 + deltaT * (i - 1), u[i - 1]) + f(t0 + deltaT * i, u[i - 1] + deltaT * f(t0 + deltaT * (i - 1), u[i - 1]));
		u[i] = u[i - 1] + 0.5 * deltaT * du;
		i++;
	}

	return u;
}

float* ODE1::modifiedEuler(float T, float deltaT)
{
	float kTot = T / deltaT, du;
	float* u = new float[kTot];
	int i;

	u[0] = cauchy;

	i = 1;
	while (i < kTot) {
		du = f(t0 + 0.5 * deltaT * (i - 1), u[i - 1] + 0.5 * deltaT * f(t0 + deltaT * (i - 1), u[i - 1]));
		u[i] = u[i - 1] + deltaT * du;
		i++;
	}

	return u;
}

float* ODE1::rungeKutta(float T, float deltaT, int s, float* a, float* b, float* c)
{
	float kTot = T / deltaT, du, uk;
	float* u = new float[kTot];
	float* k = new float[s];
	int i, j, n;

	u[0] = cauchy;

	n = 1;
	while (n < kTot) {
		i = 0;
		du = 0;
		uk = u[n - 1];
		while (i < s) {
			j = 0;
			while (j < i) {
				uk = uk + deltaT * c[i + s*j] * k[j];
				j++;
			}
			k[i] = f(t0 + deltaT * (n - 1) + b[i] * deltaT, uk);
			du = du + a[i] * k[i];
			i++;
		}
		u[n] = u[n - 1] + deltaT * du;
		n++;
	}
	
	return u;
}

