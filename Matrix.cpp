#include "Matrix.h"
#include <iostream>

void Matrix::identity()
{
	int counter = 0;
	Vector dummy = Vector(nCols);

	while (counter < nRows)
	{
		dummy.base(counter + 1);
		setColumn(counter, dummy);
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

Vector Matrix::extractColumn(int j) 
{
	Vector result = Vector(nRows);
	int i = 0;

	while (i < nRows) {
		result.setElement(i, getElement(i, j));
		i++;
	}

	return result;
}

void Matrix::setColumn(int j, Vector& column) 
{
	int i = 0;
	while (i < nRows) {
		setElement(i, j, column.getElement(i));
		i++;
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
	int j = 0;
	float factor;
	int k;

	while (j < nRows - 1) {
		int indexMaxColumn = findMaxCol(j, data, j, nRows);
		indexMaxColumn == j ? counter = counter : counter++;
		switchRows(j, indexMaxColumn);
		permutation.switchRows(j, indexMaxColumn);
		j++;
	}

	factors.push_back(counter);

	if (data[0] == 0) {
		std::cerr << "Singular matrix" << std::endl;
		exit(-1);
	}

	j = 0;
	while (j < nCols) {
		k = j + 1;
		while ( k < nRows) {
			factor = -getElement(k, j) / getElement(j, j);
			sumRows(k, j, factor);
			factors.push_back(-factor);
			k++;
		}
		j++;
	}
	return factors;
}

void Matrix::gaussianElimination()
{
	int j = 0, indexMaxColumn, k;
	float factor;

	while (j < nRows - 1) {
		indexMaxColumn = findMaxCol(j, data, j, nRows);
		switchRows(j, indexMaxColumn);
		j++;
	}

	if (data[0] == 0) {
		std::cerr << "Singular matrix" << std::endl;
		exit(-1);
	}

	j = 0;
	while (j < nCols) {
		k = j + 1;
		while (k < nRows) {
			factor = -getElement(k, j) / getElement(j, j);
			sumRows(k, j, factor);
			k++;
		}
		j++;
	}
}

Matrix Matrix::palu(Matrix& permutation, Matrix& lower)
{
	Matrix upper = Matrix(nRows);
	int j = 0, k;
	int counter = 0;
	permutation.identity();

	while (j < nRows * nCols) {
		upper.setElement(j, getElement(j));
		j++;
	}
	j = 0;

	// upper is populated via gaussian elimination
	std::vector<float> factors = upper.gaussianElimination(permutation);

	// now lower is populated 
	while (j < nRows) {
		k = j + 1;
		while (k < nRows) {
			lower.setElement(k, j, factors[counter + 1]);
			counter++;
			k++;
		}
		j++;
	}
	counter = 0;
	while (counter < nRows)
	{
		lower.setElement(counter * nRows + counter, 1);
		counter++;
	}
	return upper;
}

void Matrix::qr(Matrix& Q, Matrix& R)
{
	Vector a_0 = extractColumn(0);
	float r_00 = a_0.norm(), r_kj, r_jj;
	int j = 1, k;

	a_0 = a_0 * (1.0 / r_00);
	Q.setColumn(0, a_0);
	R.setElement(0, 0, r_00);

	Vector a_j, Q_k;

	while (j < nCols) {
		a_j = extractColumn(j);
		k = j - 1;
		while (k > -1) {
			Q_k = Q.extractColumn(k);
			r_kj = Q_k * a_j;
			R.setElement(k, j, r_kj);
			Q_k = Q_k * r_kj;
			a_j = a_j - Q_k;
			k--;
		}
		r_jj = a_j.norm();
		a_j = a_j * (1.0 / r_jj);
		Q.setColumn(j, a_j);
		R.setElement(j, j, r_jj);
		j++;
	}

}

Vector Matrix::backwardSubstitution(Vector& b)
{
	Vector x = Vector(nRows);
	float value;
	int i , j;

	x.setElement(nRows, b.getElement(nRows) / getElement(nRows, nCols));

	i = nRows - 1;
	while (i > - 1) {
		value = b.getElement(i);
		j = i + 1;
		while (j < nRows) {
			value = value - getElement(i, j) * x.getElement(j);
			j++;
		}
		x.setElement(i, value / getElement(i, i));
		i--;
	}
	return x;
}

Vector Matrix::forwardSubstitution(Vector& b)
{
	Vector x = Vector(nRows);
	float value;
	int i, j;
	
	x.setElement(0, b.getElement(0) / getElement(0));

	i = 1;
	while (i < nRows) {
		value = b.getElement(i);
		j = 0;
		while (j < i) {
			value = value - getElement(i, j) * x.getElement(j);
			j++;
		}
		x.setElement(i, value / getElement(i, i));
		i++;
	}
	return x;
}

Vector Matrix::qrSolve(Vector& b)
{
	Vector x = Vector(nRows);
	Matrix Q = Matrix(nRows), R = Matrix(nRows);

	qr(Q, R);

	Q = Q.transpose();
	b = Q * b;

	return R.backwardSubstitution(b);
}

Vector Matrix::paluSolve(Vector& b)
{
	Matrix lower = Matrix(nRows), permutation = Matrix(nRows);

	Matrix upper = palu(permutation, lower);

	b = permutation * b;

	Vector y = lower.forwardSubstitution(b);

	return upper.backwardSubstitution(y);
}

Vector Matrix::diagSolve(Vector& b) {
	Vector x = Vector(nRows);
	int i = 0;
	while (i < nRows) {
		x.setElement(i, b.getElement(i) / getElement(i, i));
		i++;
	}
	return x;
}

std::unordered_map<float, Vector> Matrix::qrEigen(int maxIterations, float tol, float* eigenValues)
{
	std::unordered_map<float, Vector> eigen;
	Matrix Q = Matrix(nRows), R = Matrix(nRows), newA = Matrix(nRows), lambda_i = Matrix(nRows);
	int numIter, i, j;
	float eigenValue, eigNorm, err = 1;
	float* eigenValues_k_1 = new float [nRows] {0.0};

	qr(Q, R);

	numIter = 0;
	while (numIter < maxIterations && err > tol) {
		newA = R * Q;
		newA.qr(Q, R);
		i = 0;
		err = 0;
		while (i < nRows) {
			eigenValues[i] = newA.getElement(i, i);
			err = err + sqrt((eigenValues[i] - eigenValues_k_1[i]) * (eigenValues[i] - eigenValues_k_1[i]));
			eigenValues_k_1[i] = eigenValues[i];
			i++;
		}
		err = err / nRows;
		numIter++;
	}

	i = 0;

	while (i < nRows) {
		eigenValue = eigenValues[i];
		lambda_i.identity();
		lambda_i = lambda_i * eigenValue;

		Matrix dummy = subtraction(lambda_i);
		dummy.gaussianElimination();
		Vector eigenVector = Vector(nRows);
		j = 0;
		while (j < nRows) {
			eigenVector.setElement(j, dummy.getElement(0, j));
			j++;
		}
		eigNorm = eigenVector.norm();
		eigenVector = eigenVector * (1.0 / eigNorm);

		eigen[eigenValue] = eigenVector;
		i++;
	}
	return eigen;
}


Vector Matrix::jacobi(Vector& b,Vector& x, float tol, int maxIter)
{
	Matrix E = Matrix(nRows), D = Matrix(nRows);
	Vector bIter = Vector(nRows), deltaB = Vector(nRows), newB = Vector(nRows);
	int i;
	float err = 1;
	float bNorm = 1.0/b.norm();

	i = 0;
	while (i < nRows) {
		D.setElement(i * nRows + i, getElement(i * nRows + i));
		i++;
	}

	E = subtraction(D);

	i = 0;
	while (err > tol && i < maxIter) {
		bIter = E * x;
		bIter = b - bIter;
		x = D.diagSolve(bIter);
		newB = product(x);
		deltaB = (b - newB);
		err = deltaB.norm() * bNorm;
		//std::cout << "Iter. n: " << i << " err: " << err << "\n";
		i++;
	}
	return x;
}

Vector Matrix::jacobi(Vector& b)
{
	int maxIter = 1e4;
	float tol = 1e-7;
	Vector x = Vector(nRows);
	Matrix E = Matrix(nRows), D = Matrix(nRows);
	Vector bIter = Vector(nRows), deltaB = Vector(nRows), newB = Vector(nRows);
	int i;
	float err = 1;
	float bNorm = 1.0 / b.norm();

	i = 0;
	while (i < nRows) {
		D.setElement(i * nRows + i, getElement(i * nRows + i));
		i++;
	}

	E = subtraction(D);

	i = 0;
	while (err > tol && i < maxIter) {
		bIter = E * x;
		bIter = b - bIter;
		x = D.diagSolve(bIter);
		newB = product(x);
		deltaB = (b - newB);
		err = deltaB.norm() * bNorm;
		std::cout << "Iter. n: " << i << " err: " << err << "\n";
		i++;
	}
	return x;
}

Vector Matrix::gaussSeidl(Vector& b, Vector& x, float tol, int maxIter)
{
	Matrix E = Matrix(nRows), F = Matrix(nRows);
	Vector bIter = Vector(nRows), deltaB = Vector(nRows), newB = Vector(nRows);
	int i, j;
	float err = 1;
	float bNorm = 1.0 / b.norm();

	i = 0;
	while (i < nRows) {
		j = 0;
		while (j < i + 1) {
			E.setElement(i, j, getElement(i, j));
			j++;
		}
		while (j < nRows) {
			F.setElement(i, j, getElement(i, j));
			j++;
		}
		i++;
	}

	i = 0;
	while (err > tol && i < maxIter) {
		bIter = F * x;
		bIter = b - bIter;
		x = E.paluSolve(bIter);
		newB = product(x);
		deltaB = (b - newB);
		err = deltaB.norm() * bNorm;
		//std::cout << "Iter. n: " << i << " err: " << err << "\n";
		i++;
	}
	return x;
}

Vector Matrix::gaussSeidl(Vector& b)
{

	int maxIter = 1e4;
	float tol = 1e-7;
	Vector x = Vector(nRows);

	Matrix E = Matrix(nRows), F = Matrix(nRows);
	Vector bIter = Vector(nRows), deltaB = Vector(nRows), newB = Vector(nRows);
	int i, j;
	float err = 1;
	float bNorm = 1.0 / b.norm();

	i = 0;
	while (i < nRows) {
		j = 0;
		while (j < i + 1) {
			E.setElement(i, j, getElement(i, j));
			j++;
		}
		while (j < nRows) {
			F.setElement(i, j, getElement(i, j));
			j++;
		}
		i++;
	}

	i = 0;
	while (err > tol && i < maxIter) {
		bIter = F * x;
		bIter = b - bIter;
		x = E.paluSolve(bIter);
		newB = product(x);
		deltaB = (b - newB);
		err = deltaB.norm() * bNorm;
		//std::cout << "Iter. n: " << i << " err: " << err << "\n";
		i++;
	}
	return x;
}



// something's wrong... i can feel it 
Vector Matrix::gradient(Vector& b, Vector& x, float tol, int maxIter)
{
	Vector z = Vector(nRows), r = Vector(nRows);
	Vector dummy = product(x);
	float normB = 1.0 / b.norm();
	int i;
	float err, alpha;

	r = b - dummy;
	err = r.norm() * normB;
	
	i = 0;
	while (i < maxIter && err > tol) {
		z = product(r);
		alpha = r * r / (r * z);
		dummy = r * alpha;
		x = x + dummy;
		dummy = z * alpha;
		z = z - dummy;

		err = r.norm() * normB;
		std::cout << "Iter. n: " << i << " err: " << err << "\n";
		i++;
	}

	return x;
}