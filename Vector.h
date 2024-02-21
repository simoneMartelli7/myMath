#pragma once
#include <iostream>

class Vector
{
private:
	int n;
	float* data;
public:
	// default constructor
	Vector() {
		this->n = 3;
		this->data = new float[n] {0};
	}
	//constructors
	Vector(int n) {
		this->n = n;
		this->data = new float[n] {0};
	};
	Vector(int n, float* input) {
		this->n = n;
		this->data = new float [n] {0};
		int j = 0;
		while (j < n) {
			data[j] = input[j];
			j++;
		}
	};
	//construct a base vector in direction i 
	//e.g. base(1) constructs [1; 0; 0] in 3d space
	void base(int i) {
		int j = 0;
		while (j < n) {
			data[j] = 0;
			j++;
		}
		data[i - 1] = 1;
	};
	Vector base(int n, int i) {
		Vector test = Vector(n);
		test.base(i);
		return test;
	};
	//print
	void print() {
		int j = 0;
		while (j < n) {
			std::cout << "| " << data[j] << " |\n";
			j++;
		}
		std::cout << "\n";
	};


	//getter
	int getN() { return this->n; };
	float* getData() { return this->data; };
	float* copyData() {
		float* result = new float[n];
		int j = 0;
		while (j < n) {
			result[j] = data[j];
			j++;
		}
		return result;
	};
	float getElement(int i) { return data[i]; };

	void setElement(int i, float value) { data[i] = value; };


	//basic operations
	Vector addition(Vector& v);
	Vector subtraction(Vector& v);
	Vector scalarProduct(float alpha);
	float dotProduct(Vector& v);


	//operators
	Vector operator+(Vector& v)   { return addition(v); };
	Vector operator-(Vector& v)   { return subtraction(v); };
	Vector operator*(float alpha) { return scalarProduct(alpha); };
	float operator*(Vector& v)    { return dotProduct(v); };
	float operator[](int i)       { return data[i]; };


	//norms
	float norm();
	float norm1();
	float normInf();
};

