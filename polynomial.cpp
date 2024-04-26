#include "polynomial.h"
#include <cmath>
#include <iostream>

float polynomial::eval(float x0) 
{
	float currentEval;
	float result = 0;
	int i = 0;

	while (i < degree+1) {
		currentEval = coeffs[i] * pow(x0, i);
		result += currentEval;
		i++;
	}

	return result;
}

void polynomial::print()
{
	int i = degree;
	std::cout << "p(x): ";
	while (i > -1) {
		std::cout << "+" << coeffs[i] << "x^" << i << " ";
		i--;
	}
	std::cout << "\n";
}

polynomial polynomial::sum(polynomial& r)
{
	int degreeMax, degreeMin, i, flag;

	if (degree > r.getDegree()) {
		degreeMax = degree;
		degreeMin = r.getDegree();
		flag = 0;
	}
	else {
		degreeMax = r.getDegree();
		degreeMin = degree;
		flag = 1;
	}

	polynomial result = polynomial(degreeMax);

	i = 0;
	while (i < degreeMin + 1) {
		result.setCoeff(i, getCoeff(i) + r.getCoeff(i));
		i++;
	}
	

	if (flag == 0) {
		while (i < degreeMax + 1) {
			result.setCoeff(i, getCoeff(i));
			i++;
		}
	}
	else {
		while (i < degreeMax + 1) {
			result.setCoeff(i, r.getCoeff(i));
			i++;
		}
	}

	return result;

}

polynomial polynomial::subtraction(polynomial& r)
{
	int degreeMax, degreeMin, i, flag;

	if (degree > r.getDegree()) {
		degreeMax = degree;
		degreeMin = r.getDegree();
		flag = 0;
	}
	else {
		degreeMax = r.getDegree();
		degreeMin = degree;
		flag = 1;
	}

	polynomial result = polynomial(degreeMax);

	i = 0;
	while (i < degreeMin + 1) {
		result.setCoeff(i, getCoeff(i) - r.getCoeff(i));
		i++;
	}


	if (flag == 0) {
		while (i < degreeMax + 1) {
			result.setCoeff(i, getCoeff(i));
			i++;
		}
	}
	else {
		while (i < degreeMax + 1) {
			result.setCoeff(i, -r.getCoeff(i));
			i++;
		}
	}

	return result;
}

polynomial polynomial::product(polynomial& r)
{
	polynomial result = polynomial(degree + r.getDegree());
	polynomial dummy = polynomial(degree + r.getDegree());
	int i, n ;
	
	i = 0;
	while (i < degree+1) {
		n = 0;
		dummy.clear();
		while (n < r.getDegree()+1) {
			dummy.setCoeff(i + n, getCoeff(i) * r.getCoeff(n));
			n++;
		}
		result = result.sum(dummy);
		i++;
	}

	return result;
}

polynomial polynomial::scalarProduct(float alpha)
{
	polynomial result = polynomial(degree);
	int i = 0;
	while (i < degree + 1) {
		result.setCoeff(i, alpha * getCoeff(i));
		i++;
	}
	return result;
}

polynomial polynomial::division(polynomial& r)
{
	if (r.getDegree() != 1) {
		std::cerr << "currently implemented for degree 1 polynomials as divisor";
		exit(-1);
	}

	polynomial result = polynomial(degree - 1);
	float k = -r.getCoeff(0), dummy;
	int i;

	i = degree-1;
	result.setCoeff(degree - 1, 1);
	dummy = 1;
	while (i > 0) {
		dummy = k * dummy;
		dummy += getCoeff(i);
		result.setCoeff(i - 1, dummy);
		i--;
	}

	if (k * dummy != -getCoeff(0)) {
		std::cerr << "the polynomial is not an exact divisor";
		exit(-1);
	}
	return result;
}
