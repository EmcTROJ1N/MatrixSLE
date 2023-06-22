#include <iostream>
#include <fstream>
#include "FileOperations.h"
#include "GaussOperations.h"

using namespace std;

unsigned i, j;
unsigned eqCount = 5;

int rows;
int cols;

double infA = -1.0;
double supA = 1.0;
double infB = 0.0;
double supB = 20.0;

double k;

double** extendedMatrix = nullptr;
double* solution = nullptr;
double* rightPart = nullptr;
double* relativeErr = nullptr;

const char* extMatrixFile = "Extended SLE matrix.txt";
const char* solutionFile = "SLE solution.txt";
const char* resultTime = "Execution time.txt";


double calculateFuncTime(double (*function)(double**, int, int), double** extendedMatrix, int rows, int cols)
{
	double functionTimeCalculate, timeRes;
	clock_t Q = clock();
	functionTimeCalculate = function(extendedMatrix, rows, cols);
	timeRes = double(clock() - Q) / CLOCKS_PER_SEC;
	return timeRes;
}

int main()
{
	system("color 70");

	clock_t Q;
	double time, resultDet, functionTimeCalculate;

	cout << "Enter matrix size: ";
	cin >> eqCount;
	rows = eqCount;
	cols = eqCount + 1;
	cout << endl;

	// Memory allocation for arrays
	rightPart = new double[rows];
	relativeErr = new double[rows];
	solution = new double[rows];
	extendedMatrix = new double* [rows];

	for (i = 0; i < rows; i++)
		extendedMatrix[i] = new double[cols];

	// Filling matrix
	k = (supA - infA) / RAND_MAX;

	for (i = 0; i < eqCount; ++i)
		for (j = 0; j < eqCount; ++j)
			extendedMatrix[i][j] = infA + k * rand();

	k = (supB - infB) / RAND_MAX;

	for (i = 0; i < rows; ++i)
		extendedMatrix[i][rows] = infB + k * rand();

	// Save matrix to the file
	if (saveMatrix(extMatrixFile, extendedMatrix, rows, cols))
		cout << "Matrix successfully saved to a file " << extMatrixFile << endl;
	else
		cout << "Something went wrong: failed to open file" << endl;

	switch (calculateGauss(extendedMatrix, rows, cols, solution))
	{
		case Unique:
		{
			if (saveSolution(solutionFile, solution, rows, cols))
				cout << "Solution successfully saved to a file " << extMatrixFile << endl;
			else
				cout << "Something went wrong: failed to open file" << endl;

			cout << "The system has a single solution!" << endl;
			cout << "The found solution is saved to a file: " << solutionFile << endl;
			cout << "Starting the test..." << endl;

			cout.setf(ios::scientific);

			cout << "    №" << '\t' << "r. column" << '\t' << "sub. r. part" << '\t' << "rel. err" << endl;

			for (i = 0; i < rows; ++i)
			{
				for (rightPart[i] = 0.0, j = 0; j < rows; ++j)
					rightPart[i] += extendedMatrix[i][j] * solution[j];

				relativeErr[i] = 1e-16 + fabs((extendedMatrix[i][rows] - rightPart[i]) / extendedMatrix[i][rows]);

				cout.width(5);
				cout << i + 1 << '\t';
				cout.precision(8);
				cout.width(8 + 5);
				cout << extendedMatrix[i][rows] << '\t' << rightPart[i] << '\t';
				cout.precision(2);
				cout.width(2 + 5);
				cout << 100 * relativeErr[i] << " %" << endl;
			}
			break;
		}

		case Empty: cout << "The system has no solutions!" << endl; break;
		case Undef: cout << "The system is undefined!" << endl; break;
		default: cout << "Something went wrong..." << endl;
	}
	cout << endl;

	

	// Calculating determinant matrix

	resultDet = calculateDeterminant(extendedMatrix, rows, cols);
	cout << "Matrix determinant: " << resultDet << endl;
	ofstream resFile(resultTime);

	resFile << "Matrix size\tCalculating time (secs)\n";

	double resTime = calculateFuncTime(calculateDeterminant, extendedMatrix, rows, cols);
	cout << fixed << setprecision(6) << "Executing time: " << resTime << endl;

	// Memory clear

	for (i = 0; i < eqCount; ++i)
		delete[] extendedMatrix[i];

	delete[] extendedMatrix;
	delete[] solution;
	delete[] rightPart;
	delete[] relativeErr;

	return 0;
}