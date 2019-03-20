
#include "pch.h" 
#include "matrix.h"
#include "stdio.h"
#include <string>
#include <algorithm>
#include <math.h>
#include "conio.h"
#include <fstream>
#include <iostream>

#pragma region someUtility
std::vector<double> load_vector(int size, std::string filename, std::vector<double> &vect) {
	FILE *in;
	fopen_s(&in, filename.c_str(), "r");
	for (int i = 0; i < size; i++) fscanf_s(in, "%lf ", &vect[i]);
	fclose(in);

	return vect;
}

void out_vector(std::vector<double> &x) {
	FILE *in;
	fopen_s(&in, "x.txt", "w");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(in, "lf", x[i]);
	}
	fclose(in);
}

int Matrix::get_dim() { return n; }

/*
double Matrix::calc_relative_discrepancy(std::vector<double> &x, std::vector<double> &f, std::vector<double> &ax) {
	ax = multiplicate_with_vector(x, ax);	    	      // Calculate Ax.
	for (int i = 0; i < ax.size(); i++) ax[i] -= f[i];    // Calculate ( Ax - f )
	double relative_discr = calc_norm(ax) / calc_norm(f); // Calculate discrepancy

	return relative_discr;
}
*/
#pragma endregion

#pragma region LUdecomposeAndSolver
void Matrix::LUdecompose() {
	// L: Lij, Lii;
	// U: Uij, 1;

/* to debug init matrix:
	mid_diag0 = std::vector<double>{ 7,1,1,1 };
	mid_diag1 = std::vector<double>{ 2,3,4,5,6 };
	mid_diag2 = std::vector<double>{ 10, 20, 30, 40, 50, 60};
	mid_diag3 = std::vector<double>{ 1, 1, 3, 2, 3};
	mid_diag4 = std::vector<double>{ 3, 2, 4, 4};
*/

	//
	mid_diag3[0] /= mid_diag2[0];
	mid_diag2[1] -= mid_diag1[0] * mid_diag3[0];

	for (int i = 1; i < n - 1; i++) {
	// caclulate elems Lij
	  //mid_diag0[i] don't change!
		mid_diag1[i] -= mid_diag0[i - 1] * mid_diag3[i - 1];

	// calculate elems Uij
		mid_diag4[i - 1] /= mid_diag2[i - 1];

		mid_diag3[i] -= mid_diag1[i-1] * mid_diag4[i-1];
		mid_diag3[i] /= mid_diag2[i];
		
	// calcuale elems Lii
		mid_diag2[i + 1] -= mid_diag0[i - 1] * mid_diag4[i - 1] + mid_diag1[i] * mid_diag3[i];

		// CHECK Lii not equal no ZERO!
	}
}

void Matrix::solveSystem(const std::vector<double>& f, std::vector<double>& q, std::vector<double>& temp) {
	// LUq = f
	
	// Uq = temp, L*temp = f, calculate temp: (forward)
	temp[0] = f[0] / mid_diag2[0];
	temp[1] = (f[1] - mid_diag1[0] * temp[0]) / mid_diag2[1];

	for (int i = 2; i < n; i++) {
		temp[i] = (f[i] - mid_diag0[i - 2] * temp[i - 2] - mid_diag1[i - 1] * temp[i - 1]) / mid_diag2[i];
	}

	// Uq = temp, calculte q: (backward)
	q[n - 1] = temp[n - 1];
	q[n - 2] = temp[n - 2] - mid_diag3[n - 2] * q[n - 1];

	for (int i = n - 3; i >= 0; i--) {
		q[i] = temp[i] - mid_diag3[i] * q[i + 1] - mid_diag4[i] * q[i + 2];
	}
}
#pragma endregion

#pragma region MultiplicateWithVector
double Matrix::calc_sum(int row, std::vector<double> &x) {
	double sum = 0;

	// Throw main diags
	if (row > 1) {
		//use 0, 1, 2 main diags
		sum += mid_diag0[row - 2] * x[row - 2];
		sum += mid_diag1[row - 1] * x[row - 1];
		sum += mid_diag2[row] * x[row];

		if (row < n - 2) { // use 3,4 diags
			sum += mid_diag3[row] * x[row + 1];
			sum += mid_diag4[row] * x[row + 2];
		}
		else if (row < n - 1) { // use all main diags
			sum += mid_diag3[row] * x[row + 1];
		}
	}
	else {
		sum += mid_diag2[row] * x[row];
		sum += mid_diag3[row] * x[row + 1];
		sum += mid_diag4[row] * x[row + 2];

		if (row == 1) {
			sum += mid_diag1[row - 1] * x[row - 1];
		}
	}

	return sum;
}

std::vector<double> Matrix::multiplicate_with_vector(std::vector<double> &x, std::vector<double> &f) {
	for (int i = 0; i < n; i++) f[i] = calc_sum(i, x);
	return f;
}
#pragma endregion

void Matrix::init(int size) {
	n = size;
	block_size = 1; // what do with block_size?
	
	// read info about accuracy of decision
	std::fstream fin(R"(input\accuracy.txt)");
	fin >> max_iter >> accuracy;
	fin.close();

	// memory allocation for Matrix
	mid_diag0 = std::vector<double>(n - 2);
	mid_diag1 = std::vector<double>(n - 1);
	mid_diag2 = std::vector<double>(n);
	mid_diag3 = std::vector<double>(n - 1);
	mid_diag4 = std::vector<double>(n - 2);
}

void Matrix::reset() {
	int i = 0;
	for (int; i < mid_diag0.size(); i++) {
		mid_diag0[i] = 0;
		mid_diag1[i] = 0;
		mid_diag2[i] = 0;
		mid_diag3[i] = 0;
		mid_diag4[i] = 0;
	}

	mid_diag1[i] = 0;
	mid_diag2[i] = 0;
	mid_diag3[i] = 0;

	mid_diag2[++i] = 0;
}

void Matrix::addLocalMatrix(const int num, const LocalMatrix& R, const LocalMatrix& M, const double dt){

	const int place = 2 * num;		// 2 == block size (0, 1, 2)

	// add elems from local matrixs to Global Matrix
	mid_diag2[place] += M[0][0] / dt + R[0][0];
	mid_diag3[place] += M[0][1] / dt + R[0][1];
	mid_diag4[place] += M[0][2] / dt + R[0][2];

	mid_diag1[place]	 += M[1][0] / dt + R[1][0];
	mid_diag2[place + 1] += M[1][1] / dt + R[1][1];
	mid_diag3[place + 1] += M[1][2] / dt + R[1][2];

	mid_diag0[place]	 += M[2][0] / dt + R[2][0];
	mid_diag1[place + 1] += M[2][1] / dt + R[2][1];
	mid_diag2[place + 2] += M[2][2] / dt + R[2][2];
}

void Matrix::setFirstBoundaryConditionsLeft() {
	mid_diag2[0] = 1;
	mid_diag3[0] = 0;
	mid_diag4[0] = 0;
}

void Matrix::setFirstBoundaryConditionsRight() {
	mid_diag0[mid_diag0.size() - 1] = 0;
	mid_diag1[mid_diag1.size() - 1] = 0;
	mid_diag2[mid_diag2.size() - 1] = 1;
}