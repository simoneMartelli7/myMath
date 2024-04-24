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
	//identity matrix
	void identity();
	//initializes the matrix with random numbers in the range [0, 1]
	void random();

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

	//BASIC OPERATIONS

	Matrix scalarProduct(float alpha);
	Matrix product(Matrix& B);
	Vector product(Vector& u);
	Vector leftProduct(Vector& u);
	Matrix addition(Matrix& B);
	Matrix subtraction(Matrix& B);
	Matrix transpose();
	Matrix diag(float value);
	Matrix diag(std::unordered_map<int, float> values);


	//OPERATORS 
	Matrix operator*(float alpha) { return scalarProduct(alpha); };
	Matrix operator*(Matrix& B)   { return product(B); };
	Vector operator*(Vector& u)   { return product(u); };
	Matrix operator+(Matrix& B)   { return addition(B); };
	Matrix operator-(Matrix& B)   { return subtraction(B); };
	Vector operator/(Vector& b)   { return solve(b); };
	

	//BASIC HANDLING
	
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
	//removes the i-th row
	Matrix removeRow(int i);
	//removes the j-th column
	Matrix removeCol(int j);
	// constructs a 1xn matrix as the transpose of a n-vector 
	//Matrix vectorTranspose(Vector& u);
	//rotates the matrix' column space of an angle theta, in 3D defaults to rotation around x-axis
	Matrix rotate(float theta);
	//rotates a matrix' 3-dimensional column space, notations is of standard euler angles
	Matrix rotate(float alpha, float beta, float gamma);
	


	// ASSESSES VARIOUS QUALITIES OF THE MATRIX 
	// 
	//returns the index of the first empty row if the matrix in not full rank
	int emptyRow();
	// returns the index of all empty rows
	int* emptyRows();
	//spectral radius of the matrix, uses qr algorithm to find eigenvalues
	float rho();
	//spectral radius, eigenvalues are passed as arguments
	float rho(float* eigenValues);
	//returns 1 if the matrix is sparse, 0 otherwise
	int isSparse();
	// retruns 1 for weak inequality, 2 for strict and 0 otherwise
	int isDiagonallyDominantRows();
	// retruns 1 for weak inequality, 2 for strict and 0 otherwise
	int isDiagonallyDominantCols();
	//returns 1 if the matrix is symmetric
	int isSymmetric();
	//returns 1 if the matrix is symmetric semi definite-positive 
	int isSymmetricDefinitePositive();


	// FACTORIZATIONS
	
	// gaussian elimination, works directly on the matrix, returns a vector containing the multiplicative factors used
	// row permutations are also done in order to have a more robus reduction, the permutation matrix doesn't need to be initialised 
	std::vector<float> gaussianElimination(Matrix& permutation);
	// gaussian elimination, works directly on the matrix 
	void gaussianElimination();
	// gaussian elimination, works directly on the matrix, doesn't raise error for singular matrix
	void gaussianSingular();
	// PA = LU factorization, returns U; P and L need to be declared but not necessarily initialised 
	Matrix palu(Matrix& permutation, Matrix& lower);
	// qr factorization, A = Q*R, both matrices should be passed as arguments, declared but not initialised 
	// uses the Gram-Schmidt alghorithm
	void qr(Matrix& Q, Matrix& R);
	//cholesky factorization, if the matrix is not symmetric, definite-positive returns an error, uses
	//the Cholesky-Crout algorithm
	Matrix cholesky();
	

	// LINEAR SYSTEMS, DIRECT 
	Vector paluSolve(Vector& b);
	Vector qrSolve(Vector& b);
	Vector backwardSubstitution(Vector& b);
	Vector forwardSubstitution(Vector& b);
	//only for diagonal matrices, if the matrix is not diagonal, the code will execute as 
	//normal but the result will be obviously wrong 
	Vector diagSolve(Vector& b);
	Vector choleskySolve(Vector& b);

	//LINEAR SYSTEMS, ITERATIVE
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
	//defaults the initial guess to a random vector, tol to 1e-7, maxIter to 1e4
	Vector gradient(Vector& b);
	// general minimal residual method with Arnoldi iteration
	//Vector gmres(Vector& b, Vector& x0, float tol, int maxIter);



	Vector solve(Vector& b);
	// EIGENVALUES 

	/*std::unordered_map<float, Vector>*/ void qrEigen(int maxIterations, float tol, float* eigenValues);


};

