#pragma once
#include <vector>
#include "matrix.h"
#include <fstream>	

struct Node {
	double x;
};

struct FiniteElem {
	int left;
	int mid;
	int right;
};



class NonlinearTask {
	using LocalMatrix = std::vector<std::vector<double>>;
	using func = std::function<double(const double&)>;
	using func2 = std::function<double(const double&, const double&)>;
	using func3d1fun = std::function<double(const double&, const double&, const double&, const func2&)>;

public:
	//NonlinearTask();
	void init();
	void setParams();

	void setFirstBoundaryConditions(const double t);

	void solve();
	void methodOfSimpleIterations();
	void linearizeLocalMatrixsAndRighPart();
	void methodOfNewton();

	void saveResult(const int timeIter, const double t);
	void resetGlobalMatrix();
	void resetGlobalF();

private:
	void calculateGlobalMatrixAndRightPart(const double t, const double dt);

	void calculateLocalMatrixOfMass( uint32_t elemNum);
	void calculateLocalMatrixOfRigid(uint32_t elemNum);
	void calculateLocalRightPart(uint32_t elemNum, const double t, const double dt);

	void addLocalRigtPartToGlobal(uint32_t elemNum);

	// Simple Iterations Method:
	bool SimpleIterationDiscrepOut();

private:
	// result output stream
	std::ofstream fout;
	std::vector<double> qExact;  // for calculate norm of error

	// time grid:
	std::vector<double> times;

	// space grid:
	std::vector<int> subareas;
	std::vector<FiniteElem> elems;
	std::vector<double> nodes;	// vector of Nodes

	// solutions:
	std::vector<double> q;
	std::vector<double> qPrev;
	std::vector<double> qPrevTime;
	std::vector<double> u;
	std::vector<double> temp;

	// local matrix of mass and rigidity, and local vector;
	LocalMatrix massLocalMatrix;
	LocalMatrix rigidLocalMatrix;
	std::vector<double> fLocal;

	// global matrix!
	std::vector<double> f;
	Matrix globalMatrix;

	// parameters of equation in different subareas:
	int amountSubareas;
	func2 uExact;
	func3d1fun fFunc;
	func2 fStart;
	std::vector<double> lambda;
	std::vector<func2> sigma;

	double leftU;		// first boundary conditions
	double rightU;

	double epsDiscrep;
};