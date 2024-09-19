#pragma once
#include <iostream>
#include <fstream>

class Vector
{
private:
	//dimension
	int n;
	//contains the values of the vector 
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
	//constructs a base vector in direction i 
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

	//fills the vector with random data in the range [0, 1]
	void random() {
		int k = 0;
		float value;
		while (k < n) {
			value = (float)(rand()) / (float)(RAND_MAX);
			setElement(k, value);
			k++;
		}
	}
	// null vector 
	void zero() 
	{
		int i = 0;
		while (i < n) {
			setElement(i, 0.0);
			i++;
		}
	}
	//print
	void print() {
		int j = 0;
		while (j < n) {
			std::cout << "| " << data[j] << " |\n";
			j++;
		}
		std::cout << "\n";
	};

	//saves to disk
	void save(std::string filename)
	{
		std::ofstream myFile;
		myFile.open(filename);

		int j = 0;
		while (j < n) {
			myFile << getElement(j) << ",\n";
			j++;
		}
		std::cout << "Vector saved in '" << filename << "'\n";
		myFile.close();
	}

	// the flag says wheter or not to overwrite the file if it already exists
	// 1 for appending, 0 otherwise 
	void save(std::string filename, int flag)
	{
		std::ofstream myFile;
		myFile.open(filename, std::ios_base::app);

		int j = 0;
		while (j < n) {
			myFile << getElement(j) << ",\n";
			j++;
		}
		std::cout << "Vector saved in '" << filename << "'\n";
		myFile.close();
	}


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

