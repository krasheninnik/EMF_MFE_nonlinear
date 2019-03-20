#pragma once
#include <vector>
#include <functional>

class Matrix {
	using LocalMatrix = std::vector<std::vector<double>>;

public:
	void init(int size);
	void reset();
	void addLocalMatrix(const int elemNum, const LocalMatrix& M, const LocalMatrix& R, const double dt);
	
	void setFirstBoundaryConditionsLeft();
	void setFirstBoundaryConditionsRight();

	int get_dim();

	double calc_sum(int, std::vector<double>&);// sum of multiplicate elems of row matrix with corresponding vector's elems
	double calc_relative_discrepancy(std::vector<double> &, std::vector<double>&, std::vector<double>&);
	std::vector<double> multiplicate_with_vector(std::vector<double>&, std::vector<double>&);

	void LUdecompose();
void solveSystem(const std::vector<double>& f, std::vector<double>& q, std::vector<double>& temp);

private: 
	static const int diags = 5;

	std::vector<double> mid_diag0;
	std::vector<double> mid_diag1;
	std::vector<double> mid_diag2;
	std::vector<double> mid_diag3;
	std::vector<double> mid_diag4;
	
	int n, block_size, max_iter;
	double accuracy;
};

