#pragma once
/*
	degree is the maximum power of the variable x
	  therefore there are n+1 coefficients ordered as increasing powers
	  e.g.: degree = 4 => ax^4 + ... + e
	  coeffs[0] = e, coeffs[degree + 1] = a
	  */
class polynomial
{
	
private:
	int degree;
	float* coeffs;
public:
	polynomial() {
		this->degree = 1;
		this->coeffs = new float[2];
		coeffs[0] = 0;
		coeffs[1] = 1;
	}
	polynomial(int degree) {
		this->degree = degree;
		this->coeffs = new float[degree + 1] {0.0};
	}
	polynomial(int degree, float* coeffsGiven) {
		this->degree = degree;
		this->coeffs = new float[degree + 1] {0.0};
		int i = 0;
		while (i < degree+1) {
			coeffs[i] = coeffsGiven[i];
			i++;
		}
	}

	//computes the polynomial in x0
	float eval(float x0);
	//printer
	void print();

	//getter
	float getDegree()     { return degree; };
	float getCoeff(int n) { return coeffs[n]; };

	//setter, only for coefficients
	void setCoeff(int n, float value) { this->coeffs[n] = value; };
	void clear()					  { this->coeffs = new float[degree + 1] {0.0}; };


	//basic operations
	polynomial sum(polynomial& r);
	polynomial subtraction(polynomial& r);
	polynomial product(polynomial& r);
	polynomial scalarProduct(float alpha);
	polynomial division(polynomial& r);

	//operators 
	polynomial operator+(polynomial& r) { return sum(r); };
	polynomial operator-(polynomial& r) { return subtraction(r); };
	polynomial operator*(polynomial& r) { return product(r); };
	polynomial operator*(float alpha) { return scalarProduct(alpha); };
	polynomial operator/(polynomial& r) { return division(r); };

};

// FIX STRANGE STUFF WITH DEGREES