#pragma once
#include <fstream>
#include <iomanip>

using namespace std;

bool saveMatrix(const char* extMatrixFile, double** extendedMatrix, int rows, int cols)
{
	ofstream file(extMatrixFile, ios::out | ios::trunc);


	if (file.is_open() == false)
		return false;
	for (unsigned i = 0; i < rows; ++i)
	{
		for (unsigned j = 0; j < cols; ++j)
			file << setw(16) << fixed << scientific << extendedMatrix[i][j];
		file << endl;
	}
	file.close();
	return true;
}

bool saveSolution(const char* solutionFile, double* solution, int rows, int cols)
{
	ofstream file(solutionFile, ios::out | ios::trunc);

	if (file.is_open() == false)
		return false;

	for (int i = 0; i < rows; ++i)
		file << std::setw(16) << fixed << scientific << solution[i] << endl;
	file.close();
}
