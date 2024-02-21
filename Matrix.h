#pragma once
#include "Vector.h"
#include <vector>
#include <unordered_map>
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
	//returns the j-th column of the matrix as a Vector 
	Vector extractColumn(int j);
	// copies the vector in the j-th column of the matrix 
	void setColumn(int j, Vector& column);


	// factorizations 
	
	// gaussian elimination, works directly on the matrix, returns a vector containing the multiplicative factors used
	// row permutations are also done in order to have a more robus reduction, the permutation matrix doesn't need to be initialised 
	std::vector<float> gaussianElimination(Matrix& permutation);
	// gaussian elimination, works directly on the matrix 
	void gaussianElimination();
	// PA = LU factorization, returns U, P and L need to be declared but not necessarily initialised 
	Matrix palu(Matrix& permutation, Matrix& lower);
	// qr factorization, A = Q*R, both matrices should be passed as arguments, declared but not initialised 
	// uses the Gram-Schmidt alghorithm
	void qr(Matrix& Q, Matrix& R);
	

	// linear systems, direct methods 
	Vector paluSolve(Vector& b);
	Vector qrSolve(Vector& b);
	Vector backwardSubstitution(Vector& b);
	Vector forwardSubstitution(Vector& b);
	//only for diagonal matrices, if the matrix is not diagonal, the code will execute as 
	//normal but the result will be obviously wrong 
	Vector diagSolve(Vector& b);

	//linear systems, iterative methods
	// jacobi method, x0 is the guess solution
	Vector jacobi(Vector& b, Vector& x0, float tol, int maxIter);
	//defaults the initial guess to the null vector, tol to 1e-7, maxIter to 1e4
	Vector jacobi(Vector& b);
	//Gauss-Seidl method, x0 is the guess solution
	Vector gaussSeidl(Vector& b, Vector& x0, float tol, int maxIter);
	//defaults the initial guess to the null vector, tol to 1e-7, maxIter to 1e4
	Vector gaussSeidl(Vector& b);
	// gradeient method
	Vector gradient(Vector& b, Vector& x0, float tol, int maxIter);


	// eigenvalues 
	std::unordered_map<float, Vector> qrEigen(int maxIterations, float tol, float* eigenValues);


};

