#include "Vector.h"
#include "Matrix.h"
#include <functional>


// FITTING

//returns the coefficient of the polynomial of degree n fitting the data y(x)
float* polynomialFitting(Vector& x, Vector& y);

//returns the value of the polynomial p(x) 
float polynomialFitting(Vector& x, Vector& y, float x0);

//returns the least square approximation to degree deg of the data 
float* leastSquareFitting(Vector& x, Vector& y, int deg);

float* lagrangePoly(Vector& x, Vector& y);


// INTEGRALS

// equally spaced nodes
Vector createNodes(float lowerLimit, float upperLimit, int nIntervals);
Vector fillNodes(Vector& x, std::function<float(float)> f);

// trapezoidal rule 
float trapz(Vector& x, Vector& y);
float trapz(float lowerLimit, float upperLimit, int N, std::function<float(float)> f);

// simpson rule, currently only implemented for equally spaced nodes 
float simpson(Vector& x, Vector& y);
float simpson(float lowerLimit, float upperLimit, int N, std::function<float(float)> f);



// NON-LINEAR EQUATIONS

// bisection method 
float bisection(float a, float b, std::function<float (float)> fx, float tol);

// Newton-Raphson method with finite differences computations 
float newtonDiscrete(float x0, std::function<float(float)> f, float tol, int maxIter);

//secant method 
float secant(float x1, float x2, std::function<float(float)> f, float tol, int maxIter);