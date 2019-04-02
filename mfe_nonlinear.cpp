#include "pch.h"
#include "mfe_nonlinear.h"
#include <fstream>
#include "matrix.h"
#include "assert.h"

void NonlinearTask::setParams() {
	lambda = std::vector<double>(amountSubareas);
	sigma = std::vector<func2>(amountSubareas);

	u = std::vector<double>(nodes.size());

	/*
		lambda[0] = 1;
		sigma[0]  = [](const double& u) {return u; };
		fFunc =		[](const double& x, const double& t, const double& u, const func& sigma) {return -2 * t + sigma(u) * x*x; };
		fStart =	[](const double& x, const double& t) {return x * x * t; };
		*/

	/* u = 5*x^4
	lambda[0] = 1;
	uExact = [](const double& x, const double& t) {return 5*x*x*x*x; };
	sigma[0] = [](const double& u, const double& x) {return 2; };
	fFunc = [](const double& x, const double& t, const double& u, const func2& sigma) {return -60*x*x; };
	*/

	lambda[0] = 1;
	uExact = [](const double& x, const double& t) {return t*t*t; };
	sigma[0] = [](const double& u, const double& x) {return 2; };
	fFunc = [](const double& x, const double& t, const double& u, const func2& sigma) {return 3*t*t*sigma(u,x); };

	//fFunc = [](const double& x, const double& t, const double& u, const func2& sigma) {return sigma(u,x) * 5 *x; };
	//fStart = [](const double& x, const double& t) {return x * t; };

	// set boundary conditions:
	leftU = uExact(nodes[0],  0);
	rightU = uExact(nodes[nodes.size() - 1], 0);
}

void NonlinearTask::init() {
#pragma region outputResultFile
	fout.open(R"(output/result.txt)");
	fout.precision(17);
#pragma endregion

#pragma region InitTimeGrid
	{
	times = std::vector<double>();

	std::fstream fin(R"(input/grid_time.txt)");
	double t0, step, coef;
	int amount, numTimeGridDividion;

	fin >> numTimeGridDividion >>  t0 >> amount >> step >> coef;

	int k = pow(2, numTimeGridDividion - 1);

	// calculate grid parameters for unevenness:
	amount *= k;
	coef = pow(coef, 1.0 / k);
	// calculate first step
	double stepsCoef = 0;
	for (int i = 0; i < k; i++) stepsCoef += pow(coef, i);	
	step /= stepsCoef;

	for (int i = 0; i < amount; i++) {
		times.push_back(t0);
		t0 += step;
		step *= coef;
	}

	fin.close();
	}
#pragma endregion

	// need change Init Grid: 1 step1 2 step1 3 
#pragma region InitSpaceGrid
	{
	nodes = std::vector<double>();
	elems = std::vector<FiniteElem>();
	subareas = std::vector<int>();

	std::fstream fin(R"(input/grid_space.txt)");

	int numSpaceGridDividion = 1;
	int numOfAreas = 1;

	fin >> numSpaceGridDividion >> numOfAreas;

	double xStart, numOfElems, step, coef;
	fin >> xStart >> numOfElems >> step >> coef;

	//const int stepsOnFiniteElem = 2;
	int k = pow(2, numSpaceGridDividion - 1);

	// calculate grid parameters for unevenness:
	numOfElems *= k;
	coef = pow(coef, 1.0 / k);
	// calculate first step
	double stepsCoef = 0;
	for (int i = 0; i < k; i++) stepsCoef += pow(coef, i);
	step /= stepsCoef;

	// there: length of 1 finite elem == 2 step
	double x = xStart;	
	nodes.push_back(x);		// add x0 in nodes
	for (int i = 0; i < numOfElems; i++) {
		nodes.push_back(x + step / 2);		// add middle node on finite elem.
		x += step;
		nodes.push_back(x);

		step *= coef;					// change step
	}

	// fill elems array:
	int j = 0;
	for (; j < nodes.size() - 2; j+=2) {
		elems.push_back(FiniteElem{ j, j + 1, j + 2});
	}

	assert(j == nodes.size() - 1);

	// fill subareas:
	int  numOfFiniteElems = 0;
	int sum = 0;

	fin >> amountSubareas;
	for (int i = 0; i < amountSubareas; i++) {
		fin >> numOfFiniteElems;
		numOfFiniteElems *= k; // consideration of grid dividion
		sum += numOfFiniteElems;
		for(int j = 0; j < numOfFiniteElems; j++) subareas.push_back(0);
	}	
	
	assert(sum == elems.size());

	// init vector of params of equals in subareas:
	lambda = std::vector<double>(amountSubareas);
	sigma = std::vector<func2>(amountSubareas);

	fin.close();
	}
#pragma endregion

#pragma region MatrixInit
	globalMatrix = Matrix();
	
	const int matrixDim = 3 + (elems.size() - 1) * 2;
	globalMatrix.init(matrixDim); 

#pragma endregion

#pragma region MemoryAllocation

	q = std::vector<double>(nodes.size());
	qPrev = std::vector<double>(nodes.size());
	qPrevTime = std::vector<double>(nodes.size());
	qExact = std::vector<double>(nodes.size());;


	// there should be vector<vvector<doub>> results -> for save q on each step!
	temp = std::vector<double>(nodes.size());

	// local matrix of mass and rigidity, and local vector;
	const int pointsToFiniteElem = 3;
	massLocalMatrix = LocalMatrix(pointsToFiniteElem);
	for (auto& el : massLocalMatrix) el = std::vector<double>(pointsToFiniteElem);
	
	rigidLocalMatrix = LocalMatrix(pointsToFiniteElem);
	for (auto& el : rigidLocalMatrix) el = std::vector<double>(pointsToFiniteElem);

	fLocal = std::vector<double>(pointsToFiniteElem);

	// global matrix!
	f = std::vector<double>(nodes.size());;
	
	//TapeMatrix globalMatirx;
#pragma endregion
}

void NonlinearTask::calculateGlobalMatrixAndRightPart(const double t, const double dt) {
	// clear global matrix and vector:
	resetGlobalMatrix();
	resetGlobalF();

	// calculate new global matrix and vector:
	for (int elemNum = 0; elemNum < elems.size(); elemNum++) {
		calculateLocalMatrixOfMass(elemNum);	// depends on q
		calculateLocalMatrixOfRigid(elemNum);
		calculateLocalRightPart(elemNum, t, dt);	// depends on q

		// if method Newton'a need linearizate local matrix
		globalMatrix.addLocalMatrix(elemNum, rigidLocalMatrix, massLocalMatrix, dt);
		addLocalRigtPartToGlobal(elemNum);
	}

	// set boundary conditions
	setFirstBoundaryConditions(t);
}
#include <iostream>
void NonlinearTask::solve() {
	// need input things:
	epsDiscrep = 1e-15;

	double t = 0;
	// calculate u0:
	for (uint32_t i = 0; i < u.size(); i++) {
		q[i] = qPrevTime[i] = uExact(nodes[i], times[0]);
	}

	// approximation by time:
	for (uint32_t i = 1; i < times.size(); i++) {
		const double dt = times[i] - times[i - 1];
		const double t = times[i];
		calculateGlobalMatrixAndRightPart(t, dt);

		// solve nonlinear system of equations: A(qi)qi = b(qi)
		// there mathod of simple itertations: A(qi-1)qi = b(qi-1)

		bool calculating = true;
		while (calculating) {
			globalMatrix.LUdecompose();
			globalMatrix.solveSystem(f, q, temp);
			
			// there we have new q.
			// calculating A(qNew) and f(qNew )
			calculateGlobalMatrixAndRightPart(t, dt);
			
			if (SimpleIterationDiscrepOut()) {
				calculating = false;
				saveResult(i, t);

				std::cout << "+" << std::endl;
			}
			else {
				// prepare next iteration:
				std::swap(q, qPrev);
			}
		}
		
		// prepare next iteration:
		//std::swap(q, qPrevTime);
		qPrevTime = q;  
	}

}

void NonlinearTask::calculateLocalMatrixOfMass(uint32_t elemNum){
	// approximate sigma U part
	const auto& elem = elems[elemNum];
	const auto k = (nodes[elem.right] - nodes[elem.left]);

						  // sigma_k(u)
	const auto& sigmaLocal = sigma[subareas[elemNum]];
	const double coef0 = k * sigmaLocal(q[elem.left], nodes[elem.left]);		// sigma(u[elem.first)]; * multiply to step. [nodes[elem.right] - nodes[elem.left]);
	const double coef1 = k * sigmaLocal(q[elem.mid], nodes[elem.mid]);		// sigma(u[elem.second)]'
	const double coef2 = k * sigmaLocal(q[elem.right], nodes[elem.right]);

	massLocalMatrix[0][0] = coef0 * 0.09285714285714286 + coef1 * 0.04761904761904761 + coef2 * -7.1428571428571415e-3;
	massLocalMatrix[0][1] = coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905;
	massLocalMatrix[0][2] = coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3;
							 
	massLocalMatrix[1][0] = coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905;
	massLocalMatrix[1][1] = coef0 * 0.0380952380952381	+ coef1 * 0.4571428571428571 + coef2 * 0.03809523809523809;
	massLocalMatrix[1][2] = coef0 * -0.01904761904761905 +coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761;
							
	massLocalMatrix[2][0] = coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3;
	massLocalMatrix[2][1] = coef0 * -0.01904761904761905	+ coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761;
	massLocalMatrix[2][2] = coef0 * -7.1428571428571415e-3 + coef1 * 0.04761904761904761 + coef2 * 0.09285714285714286;
}

void NonlinearTask::calculateLocalMatrixOfRigid(uint32_t elemNum) {
	// approximate div(lambda grad) part:								
	const auto& elem = elems[elemNum];	
	const double coef = lambda[subareas[elemNum]] / ( 3 * (nodes[elem.right] - nodes[elem.left])); 

	rigidLocalMatrix[0][0] = coef *  7;
	rigidLocalMatrix[0][1] = coef * -8;
	rigidLocalMatrix[0][2] = coef *  1;
							 
	rigidLocalMatrix[1][0] = coef * -8;
	rigidLocalMatrix[1][1] = coef * 16;
	rigidLocalMatrix[1][2] = coef * -8;
							 
	rigidLocalMatrix[2][0] = coef *  1;
	rigidLocalMatrix[2][1] = coef * -8;
	rigidLocalMatrix[2][2] = coef *  7;
}

void NonlinearTask::calculateLocalRightPart(uint32_t num, const double t, const double dt) {
	// f = f + (1/dt * M * q_j-1)
	const int size = 3;
	
	const int block = 2;
	int place = num * block;

	for (int i = 0; i < size; i++) {
		fLocal[i] = 0;
		for (int j = 0; j < size; j++) {
			fLocal[i] += massLocalMatrix[i][j] * qPrevTime[place + j];
		}
		fLocal[i] /= dt;
	}

	const auto& sigmaLocal = sigma[subareas[num]];

											// there should be q or qPrevTime !?
											// => q. because solve the equation in "q" time.
	const double f0 = fFunc(nodes[place], t, q[place], sigmaLocal);
	const double f1 = fFunc(nodes[place + 1], t, q[place + 1], sigmaLocal);
	const double f2 = fFunc(nodes[place + 2], t, q[place + 2], sigmaLocal);

	const auto& elem = elems[num];
	const auto k = (nodes[elem.right] - nodes[elem.left]) / 30;

	fLocal[0] += k * (4*f0 + 2*f1 - f2) ;
	fLocal[1] += k * (2*f0 + 16*f1 + 2*f2);
	fLocal[2] += k * (-f0 + 2*f1 + 4*f2);
}

void NonlinearTask::addLocalRigtPartToGlobal(uint32_t num) {
	const int block = 2;
	int place = num * block;

	for (int i = 0; i <= block; i++, place++) f[place] += fLocal[i];
}

void NonlinearTask::setFirstBoundaryConditions(const double t) {
	// set first bounday conditions in left side:
	globalMatrix.setFirstBoundaryConditionsLeft();
	f[0] = uExact(nodes[0], t);

	// set first bounday conditions in right side:
	globalMatrix.setFirstBoundaryConditionsRight();
	f[f.size() - 1] = uExact(nodes[nodes.size() - 1], t);
}

void vectorSubtraction(std::vector<double>& result, const std::vector<double>& a){
	for (int i = 0; i < result.size(); i++) result[i] -= a[i];
}

double calcNorm(const std::vector<double> &x) {
	double norm = 0;
	for (int i = 0; i < x.size(); i++) {
		norm += x[i] * x[i];
	}
	norm = sqrt(norm);
	return norm;
}

//#include <iostream>
bool NonlinearTask::SimpleIterationDiscrepOut() {
	// || A(qi) * qi - b(qi) || / || b(qi) || < eps => out:

	temp = globalMatrix.multiplicate_with_vector(q, temp);
	vectorSubtraction(temp, f);

	double resultNorm = calcNorm(temp) / calcNorm(f);
	//std::cout << resultNorm << std::endl;

	return resultNorm < epsDiscrep;
}

void NonlinearTask::resetGlobalMatrix() {
	globalMatrix.reset();
}

void NonlinearTask::resetGlobalF() {
	for (int i = 0; i < f.size(); i++) f[i] = 0;
}

void NonlinearTask::saveResult(const int timeIter, const double t) {
	double sum = 0;
	for (int i = 0; i < qExact.size(); i++) {
		qExact[i] = uExact(nodes[i], t);
		sum += (qExact[i] - q[i])*(qExact[i] - q[i]);
	}
	sum = sqrt(sum);

	fout << "TimeIteration: " << timeIter << "\tTime: " << t << "\tNormOfError: " << sum << std::endl;
	for (const auto& resExact : qExact) fout << resExact << " ";  fout << std::endl;
	for (const auto& res : q) fout << res << " "; fout << std::endl;
	fout << std::endl;
}