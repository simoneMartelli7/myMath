#include "ODE.h"

Vector* ODE::explicitEulerNoForcing(float deltaT, float T)
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

Vector* ODE::explicitEulerForcing(float deltaT, float T)
{
	float kTot = T / deltaT;
	Vector du, dummy;
	Vector* u = new Vector[kTot];
	int i;

	u[0] = cauchy;

	i = 1;
	while (i < kTot) {
		du = A * u[i - 1];
		dummy = g.eval(t0 + deltaT * (i - 1));
		du = du + dummy;
		du = du * deltaT;
		u[i] = u[i - 1] + du;
		i++;
	}

	return u;
}


Vector* ODE::heunNoForcing(float deltaT, float T)
{
	float kTot = T / deltaT;
	Vector du1, du2, u2;
	Vector* u = new Vector[kTot];
	int i;

	u[0] = cauchy;

	i = 1;
	while (i < kTot) {
		//f(t_k, y_k)
		du1 = A * u[i - 1];

		//f(t_k, y_k + hf(t_k, y_k))
		du2 = du1 * deltaT;
		u2 = u[i - 1] + du2;
		du2 = A * u2;

		du1 = du1 + du2;
		du1 = du1 * (0.5 * deltaT);
		u[i] = u[i - 1] + du1;

		i++;
	}

	return u;
}

Vector* ODE::heunForcing(float deltaT, float T)
{
	float kTot = T / deltaT;
	Vector du1, dummy, du2, u2;
	Vector* u = new Vector[kTot];
	int i;

	u[0] = cauchy;

	i = 1;
	while (i < kTot) {
		//f(t_k, y_k)
		Vector u_k = u[i - 1];

		float imDumb = 0;

		du1 = A * u_k;
		dummy = g.eval(t0 + deltaT * (i - 1));
		du1 = du1 + dummy;

		//f(t_k, y_k + hf(t_k, y_k))
		du2 = du1 * deltaT;
		du2 = du2 + dummy;
		u2 = u_k + du2;
		du2 = A * u2;
		dummy = g.eval(t0 + deltaT * i);
		du2 = du2 + dummy;
		
		du1 = du1 + du2;
		du1 = du1 * (0.5 * deltaT);
		u[i] = u_k + du1;

		i++;
	}

	return u;
}


