#pragma once
#include "Vector.h"
#include <vector>
class Matrix
{
private:
	int nRows;
	int nCols;
	float* data;
public:
	//default constructor
	Matrix() {
		this->nRows = 3;
		this->nCols = 3;
		this->data = new float[nRows * nCols] {0.0};
	}
	//constructor
	Matrix(int nRows, int nCols) {
		this->nRows = nRows;
		this->nCols = nCols;

		this->data = new float[nRows * nCols] {0.0};
	}
	//square matrix
	Matrix(int n) {
		this->nRows = n;
		this->nCols = n;
		this->data = new float[n * n] {0.0};
	}
	//identity
	void identity();

	//getter
	int getCols() { return nCols; };
	int getRows() { return nRows; };
	float* getData() { return data; };
	float getElement(int i, int j) { return data[i * nCols + j]; };
	float getElement(int i) { return data[i]; };

	//setter
	void setElement(int i, int j, float value) { data[i * nCols + j] = value; };
	void setElement(int i, float value) { data[i] = value; };
	//printer
	void print();

	//basic operations
	Matrix scalarProduct(float alpha);
	Matrix product(Matrix& B);
	Vector product(Vector& u);
	Matrix addition(Matrix& B);
	Matrix subtraction(Matrix& B);
	Matrix transpose();

	//operators
	Matrix operator*(float alpha) { return scalarProduct(alpha); };
	Matrix operator*(Matrix& B)   { return product(B); };
	Vector operator*(Vector& u)   { return product(u); };
	Matrix operator+(Matrix& B)   { return addition(B); };
	Matrix operator-(Matrix& B)   { return subtraction(B); };
	

	//handling operations
	
	//switches row1 with row2  
	void switchRows(int row1, int row2);
	//switches col1 with col2
	void switchCols(int col1, int col2);
	// row1 = row1+alpha*row2
	void sumRows(int row1, int row2, float alpha);


	// factorizations 
	std::vector<float> gaussianElimination(Matrix& permutation);


};

