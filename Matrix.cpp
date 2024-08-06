#include "Matrix.h"
#include <iostream>
#include <fstream>

//BASIC OPERATIONS

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

void Matrix::random() {
	int k = 0;
	float value;
	while (k < nRows*nCols) {
		value = (float)(rand()) / (float)(RAND_MAX);
		setElement(k, value);
		k++;
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

void Matrix::save(std::string filename)
{
	std::ofstream myFile;
	myFile.open(filename);

	int i, j;

	i = 0;
	while (i < nRows) {
		j = 0;
		while (j < nCols) {
			myFile << getElement(i, j) << ", ";
			j++;
		}
		myFile << "\n";
		i++;
	}
	std::cout << "Matrix saved in '" << filename << "'\n";
	myFile.close();
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

Vector Matrix::leftProduct(Vector& u)
{
	if (nRows != u.getN()) {
		std::cerr << "incompatible dimensions";
		exit(-1);
	}
	Vector result = Vector(nCols);
	int i, j;
	float v_j;

	j = 0;
	while (j < nCols) {
		v_j = 0;
		i = 0;
		while (i < nRows) {
			v_j += u.getElement(i) * getElement(i, j);
			j++;
		}
		result.setElement(i, v_j);
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
	Matrix result = Matrix(nCols, nRows);

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

Matrix Matrix::diag(float value) 
{
	identity();
	return scalarProduct(value);
}

Matrix Matrix::diag(std::unordered_map<int, float> values)
{
	int counter;
	std::unordered_map<int, float>::iterator iterator;
	Matrix result = Matrix(nRows, nCols);

	counter = 0;
	while (counter < nRows) {
		iterator = values.begin();
		while (iterator != values.end()) {
			result.setElement(counter + nCols * counter + iterator->first, iterator->second);
			iterator++;
		}
		counter++;
	}
	return result;
}



//BASIC HANDLING

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

Matrix Matrix::removeRow(int i)
{
	Matrix result = Matrix(nRows - 1, nCols);
	int k = 0;
	
	while (k < nCols * i) {
		result.setElement(k, getElement(k));
		k++;
	}

	k = k + nCols;
	while (k < nCols * nRows) {
		result.setElement(k - nCols, getElement(k));
		k++;
	}
	return result;
}

Matrix Matrix::removeCol(int j)
{
	Matrix dummy = transpose();

	return dummy.removeRow(j).transpose();
}

Matrix Matrix::rotate(float theta)
{
	if (nCols == 2)
	{
		Matrix rotation = Matrix(2);
		rotation.setElement(0,  cos(theta));
		rotation.setElement(1, -sin(theta));
		rotation.setElement(2,  sin(theta));
		rotation.setElement(3,  cos(theta));

		return product(rotation);
	}
	else if (nCols == 3)
	{
		Matrix rotation = Matrix(3);
		rotation.setElement(0, 1);
		rotation.setElement(4,  cos(theta));
		rotation.setElement(5, -sin(theta));
		rotation.setElement(7,  sin(theta));
		rotation.setElement(8,  cos(theta));

		return product(rotation);
	}

}

Matrix Matrix::rotate(float alpha, float beta, float gamma)
{

	if (nCols != 3) {
		std::cerr << "Function only implemented for 3x3 matrices";
		exit(-1);
	}

	Matrix rotation = Matrix(3), r_i = Matrix(3);

	rotation.setElement(0, cos(alpha));
	rotation.setElement(1, -sin(alpha));
	rotation.setElement(3, sin(alpha));
	rotation.setElement(4, cos(alpha));
	rotation.setElement(8, 1);

	r_i.setElement(0, cos(beta));
	r_i.setElement(2, sin(beta));
	r_i.setElement(4, 1);
	r_i.setElement(6, -sin(beta));
	r_i.setElement(8, cos(beta));

	rotation = rotation * r_i;

	r_i.setElement(0, 1);
	r_i.setElement(4, cos(gamma));
	r_i.setElement(5, -sin(gamma));
	r_i.setElement(7, sin(gamma));
	r_i.setElement(8, cos(gamma));

	rotation = rotation * r_i;

	return product(rotation);
}




// ASSESSES VARIOUS QUALITIES OF THE MATRIX 

int Matrix::emptyRow()
{
	int i, j;
	int emptyRow;

	i = 0;
	while (i < nRows) {
		j = 0;
		emptyRow = 0;
		while (j < nCols) {
			emptyRow += getElement(i, j);
			j++;
		}
		if (emptyRow == 0) {
			return i;
		}
		i++;
	}

	std::cerr << "The matrix is full rank, use paluSolve or qrSolve";
	exit(-1);
}

int* Matrix::emptyRows()
{
	int i, j;
	int emptyRow;
	int* emptyRows = new int[nRows] {0};

	i = 0;
	while (i < nRows) {
		j = 0;
		emptyRow = 0;
		while (j < nCols) {
			emptyRow += getElement(i, j);
			j++;
		}
		if (emptyRow == 0) {
			emptyRows[i] = 1;;
		}
		i++;
	}
	return emptyRows;
}

float Matrix::rho()
{
	float* eigen = new float[nRows] {0.0};
	float max;
	int i;

	qrEigen(1e3, 1e-6, eigen);

	i = 1;
	max = abs(eigen[0]);
	while (i < nRows) {
		if (max < abs(eigen[i])) {
			max = abs(eigen[i]);
		}
		i++;
	}
	return max;
}

float Matrix::rho(float* eigenValues)
{
	float max;
	int i;

	i = 1;
	max = abs(eigenValues[0]);
	while (i < nRows) {
		if (max < abs(eigenValues[i])) {
			max = abs(eigenValues[i]);
		}
		i++;
	}
	return max;
}

int Matrix::isSparse()
{
	int counter, i;
	float k = 0.2; // this may need some tweaking 

	counter = 0;
	i = 0;
	while (i < nRows) {
		if (getElement(i) != 0) {
			counter++;
		}
		i++;
	}

	return k * nRows* nCols > counter ? 1 : 0;
}

int Matrix::isDiagonallyDominantRows()
{
	float sum;
	int i, j, flag = 2;

	i = 0;
	while (i < nRows) {
		j = 0;
		sum = 0;
		while (j != i) {
			sum += abs(getElement(i, j));
			j++;
		}
		j++;
		while (j < nCols) {
			sum += abs(getElement(i, j));
			j++;
		}

		if (abs(getElement(i, i)) > sum && flag == 2) {
			i++;
		}
		else if (abs(getElement(i, i)) >= sum && (flag == 2 || flag == 1)) {
			flag = 1;
			i++;
		}
		else
			return 0;
	}
	return flag;
}

int Matrix::isDiagonallyDominantCols()
{
	float sum;
	int i, j, flag = 2;

	i = 0;
	while (i < nCols) {
		j = 0;
		sum = 0;
		while (j != i) {
			sum += abs(getElement(j, i));
			j++;
		}
		j++;
		while (j < nRows) {
			sum += abs(getElement(j, i));
			j++;
		}

		if (abs(getElement(i, i)) > sum && flag == 2) {
			i++;
		}
		else if (abs(getElement(i, i)) >= sum && (flag == 2 || flag == 1)) {
			flag = 1;
			i++;
		}
		else
			return 0;
	}
	return flag;
}

int Matrix::isSymmetric()
{
	int i, j;

	i = 0;
	while (i < nRows) {
		j = i+1;
		while (j < nCols) {
			if (abs(getElement(i, j) - getElement(j, i)) > 1e-7) {
				return 0;
			}
			j++;
		}
		i++;
	}
	return 1;
}

int Matrix::isSymmetricDefinitePositive()
{
	if (isSymmetric()) {
		float* eigen = new float[nRows] {0.0};
		qrEigen(1e3, 1e-7, eigen);
		int i = 0;
		while (i < nRows) {
			if (eigen[i] < 0) {
				//early exit when one eigenvalue is negative 
				return 0;
			}
			i++;
		}
	}
	return 1;
}




// FACTORIZATIONS

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

void Matrix::permutation(Matrix& permutation)
{
	permutation.identity();
	float* permutationData = permutation.getData();

	int j = 0;

	while (j < nRows - 1) {
		int indexMaxColumn = findMaxCol(j, data, j, nRows);
		switchRows(j, indexMaxColumn);
		permutation.switchRows(j, indexMaxColumn);
		j++;
	}
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

void Matrix::gaussianSingular()
{
	int j = 0, indexMaxColumn, k;
	float factor;

	while (j < nRows - 1) {
		indexMaxColumn = findMaxCol(j, data, j, nRows);
		switchRows(j, indexMaxColumn);
		j++;
	}

	int* zeroRows = emptyRows();
	
	
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

Matrix Matrix::cholesky()
{
	if (!isSymmetricDefinitePositive()) {
		std::cerr << "The matrix is not symmetric definite-positive";
		exit(-1);
	}

	Matrix L = Matrix(nRows);
	int i, j, k;
	float sum;

	j = 0;
	while (j < nRows) {
		sum = 0;
		k = 0;
		while (k < j) {
			sum += L.getElement(j, k) * L.getElement(j, k);
			k++;
		}
		L.setElement(j, j, sqrt(getElement(j, j) - sum));

		i = j + 1;
		while (i < nRows) {
			sum = 0;
			k = 0;
			while (k < j) {
				sum += L.getElement(i, k) * L.getElement(j, k);
				k++;
			}
			L.setElement(i, j, 1.0 / L.getElement(j, j) * (getElement(i, j) - sum));
			i++;
		}
		j++;
	}
	return L;
}



// LINEAR SYSTEMS, DIRECT 

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

Vector Matrix::choleskySolve(Vector& b)
{
	Vector x = Vector(nRows), y = Vector(nRows);
	Matrix L = cholesky(), Lt = L.transpose();

	y = L.forwardSubstitution(b);
	x = Lt.backwardSubstitution(y);

	return x;
}


//LINEAR SYSTEMS, ITERATIVE

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

Vector Matrix::gradient(Vector& b, Vector& x, float tol, int maxIter)
{
	Vector t = product(x);
	Vector r = b - t, p = r, rPlus = Vector(nRows), dummy = Vector(nRows);
	int k;
	float alpha, beta;
	float err = r.norm();

	k = 0;
	while (k < maxIter && r.norm() > tol) {
		t = product(p);
		alpha = r.dotProduct(r) / p.dotProduct(t);
		dummy = p * alpha;
		x = x + dummy;
		dummy = t * alpha;
		rPlus = r - dummy;
		beta = rPlus.dotProduct(rPlus) / (r.dotProduct(r));
		dummy = p * beta;
		p = rPlus + dummy;
		r = rPlus;
		k++;
	}
	return x;
}

Vector Matrix::gradient(Vector& b)
{
	Vector x = Vector(nRows);
	x.random();
	Vector t = product(x);
	Vector r = b - t, p = r, rPlus = Vector(nRows), dummy = Vector(nRows);
	int k, maxIter = 1e4;
	float alpha, beta, tol = 1e-7;
	float err = r.norm();

	k = 0;
	while (k < maxIter && r.norm() > tol) {
		t = product(p);
		alpha = r.dotProduct(r) / p.dotProduct(t);
		dummy = p * alpha;
		x = x + dummy;
		dummy = t * alpha;
		rPlus = r - dummy;
		beta = rPlus.dotProduct(rPlus) / (r.dotProduct(r));
		dummy = p * beta;
		p = rPlus + dummy;
		r = rPlus;
		k++;
	}
	return x;
}

/*Vector Matrix::gmres(Vector& b, Vector& x0, float tol, int maxIter)
{
	Vector r = Vector(nRows), dummy = product(x0);
	float err = 1, beta;
	float* cs = new float[maxIter];
	float* sn = new float[maxIter];

	r = b - dummy;
	err = r.norm();
	beta = err;

}*/

// NOT ACTUALLY IMPLEMENTED
Vector Matrix::solve(Vector& b) 
{
	if (isSparse()) {
		return gradient(b);
	}
	else {
		return paluSolve(b);
	}
}




// EIGENVALUES

void Matrix::qrEigen(int maxIterations, float tol, float* eigenValues)
{
	std::unordered_map<float, Vector> eigen;
	Matrix Q = Matrix(nRows), R = Matrix(nRows), newA = Matrix(nRows), eigenVectors = Matrix(nRows);
	int numIter, i;
	float err = 1;
	float* eigenValues_k_1 = new float [nRows] {0.0};

	qr(Q, R);
	eigenVectors.identity();

	numIter = -1;
	while (numIter++ < maxIterations && err > tol) {
		newA = R * Q;
		newA.qr(Q, R);
		i = -1;
		err = 0;
		//eigenVectors = eigenVectors * Q;
		while (i++ < nRows) {
			eigenValues[i] = newA.getElement(i, i);
			err = err + sqrt((eigenValues[i] - eigenValues_k_1[i]) * (eigenValues[i] - eigenValues_k_1[i]));
			eigenValues_k_1[i] = eigenValues[i];
			//i++;
		}
		err = err / nRows;
		//numIter++;
	}

	/*i = 0;
	while (i < nCols) {
		eigen[eigenValues[i]] = eigenVectors.extractColumn(i);
		i++;
	}*/
	//return eigen;
}