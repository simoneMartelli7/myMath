#include "Matrix.h"
#include <iostream>

void Matrix::identity()
{
	int counter = 0;
	while (counter < nRows)
	{
		data[counter * nRows + counter] = 1;
		counter++;
	}
}

void Matrix::print()
{
	int i = 0;
	int j = 0;
	while (i < nRows) {
		std::cout << "| ";
		j = 0;
		while (j < nCols) {
			std::cout << getElement(i, j) << " ";
			j++;
		}
		std::cout << "|\n";
		i++;
	}
	std::cout << "\n";
}

Matrix Matrix::scalarProduct(float alpha) {
	Matrix result = Matrix(nRows, nCols);
	int i = 0;
	while (i < nRows * nCols) {
		result.setElement(i, getElement(i) * alpha);
		i++;
	}
	return result;
}

Matrix Matrix::product(Matrix& B) 
{
	int m1 = nRows;
	int n1 = nCols;

	int m2 = B.getRows();
	int n2 = B.getCols();

	if (n1 != m2) {
		std::cerr << "Error: matrix multiplication is defined only in the following case (m x n) * (n x p)" << std::endl;
		exit(-1);
	}

	Matrix result = Matrix(m1, n2);

	int i = 0, j = 0, k = 0;

	while (i < m1) {
		j = 0;
		while (j < n2) {
			float c_ij = 0;
			k = 0;
			while (k < n1) {
				c_ij += getElement(i, k) * B.getElement(k, j);
				k++;
			}
			result.setElement(i, j, c_ij);
			j++;
		}
		i++;
	}
	return result;
}

Vector Matrix::product(Vector& u)
{
	if (nCols != u.getN()) {
		std::cerr << "incompatible dimensions";
		exit(-1);
	}
	Vector result = Vector(nRows);
	int i = 0, j = 0;
	float v_i = 0;
	while (i < nRows) {
		v_i = 0;
		j = 0;
		while (j < nCols) {
			v_i += getElement(i, j) * u.getElement(j);
			j++;
		}
		result.setElement(i, v_i);
		i++;
	}
	return result;
}

Matrix Matrix::addition(Matrix& B)
{
	if (nRows != B.getRows() || nCols != B.getCols()) {
		std::cerr << "incompatible dimensions\n";
		exit(-1);
	}
	Matrix result = Matrix(nRows, nCols);

	int i = 0;
	while (i < nRows * nCols) {
		result.setElement(i, getElement(i) + B.getElement(i));
		i++;
	}
	return result;
}

Matrix Matrix::subtraction(Matrix& B)
{
	if (nRows != B.getRows() || nCols != B.getCols()) {
		std::cerr << "incompatible dimensions\n";
		exit(-1);
	}
	Matrix result = Matrix(nRows, nCols);

	int i = 0;
	while (i < nRows * nCols) {
		result.setElement(i, getElement(i) - B.getElement(i));
		i++;
	}
	return result;
}

Matrix Matrix::transpose()
{
	Matrix result = Matrix(nRows, nCols);

	int i = 0, j;
	while (i < nCols) {
		j = 0;
		while (j < nRows) {
			result.setElement(i, j, getElement(j, i));
			j++;
		}
		i++;
	}
	return result;
}

void Matrix::switchRows(int row1, int row2)
{
	if (row1 >= nRows || row2 >= nRows || row1 < 0 || row2 < 0) {
		std::cerr << "Input error: selected row(s) doesn't belong to the matrix" << std::endl;
		exit(-1);
	}
	int j = 0;
	float copy1, copy2;
	while (j < nCols) {
		copy1 = getElement(row1, j);
		copy2 = getElement(row2, j);

		setElement(row1, j, copy2);
		setElement(row2, j, copy1);
		j++;
	}
}

void Matrix::switchCols(int col1, int col2)
{
	if (col1 >= nRows || col2 >= nRows || col1 < 0 || col2 < 0) {
		std::cerr << "Input error: selected row(s) doesn't belong to the matrix" << std::endl;
		exit(-1);
	}
	int j = 0;
	float copy1, copy2;
	while (j < nRows) {
		copy1 = getElement(j, col1);
		copy2 = getElement(j, col2);

		setElement(j, col1, copy2);
		setElement(j, col2, copy1);
		j++;
	}
}

void Matrix::sumRows(int row1, int row2, float alpha)
{
	int j = 0;
	while (j < nCols) {
		setElement(row1, j, getElement(row1, j) + alpha * getElement(row2, j));
		j++;
	}
}

int findMaxCol(int col, float* A, int firstRow, int n) {
	float max = A[firstRow * n + col];
	int maxIndex = firstRow;
	int i = firstRow + 1;

	while (i < n) {
		if (abs(A[firstRow * n + col]) > max) {
			max = abs(A[firstRow * n + col]);
			maxIndex = i;
		}
		i++;
	}
	return maxIndex;
}

std::vector<float> Matrix::gaussianElimination(Matrix& permutation) {

	permutation.identity();
	float* permutationData = permutation.getData();

	int counter = 0;
	std::vector<float> factors;

	for (int j = 0; j < nRows - 1; j++) {
		int indexMaxColumn = findMaxCol(j, data, j, nRows);
		indexMaxColumn == j ? counter = counter : counter++;
		switchRows(j, indexMaxColumn);
		permutation.switchRows(j, indexMaxColumn);
	}
	factors.push_back(counter);

	if (data[0] == 0) {
		std::cerr << "Singular matrix" << std::endl;
		exit(-1);
	}

	for (int i = 0; i < nCols; i++) {

		for (int k = i + 1; k < nRows; k++) {
			float factor = -getElement(k, i) / getElement(i, i);
			sumRows(k, i, factor);
			factors.push_back(-factor);
		}
	}
	return factors;
}