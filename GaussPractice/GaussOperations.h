#pragma once
#include <cmath>
#include <iomanip>

using namespace std;

enum GaussSolutionResult
{
	Unique,
	Empty,
	Undef,
	Error
};

GaussSolutionResult calculateGauss(double** extendedMatrix, int rows, int cols, double* solution)
{
	GaussSolutionResult res = Unique;

	double** tmpMatrix = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		tmpMatrix[i] = new double[cols];
		memcpy(tmpMatrix[i], extendedMatrix[i], cols * sizeof(double));
	}

	for (unsigned k = 0; k < rows - 1; k++)
	{
		double maxValue = abs(tmpMatrix[k][k]);
		unsigned maxIdx = k;

		for (int i = k + 1; i < rows; i++)
		{

			if (abs(tmpMatrix[i][k]) > maxValue)
			{
				maxValue = abs(tmpMatrix[i][k]);
				maxIdx = i;
			}
		}

		if (maxIdx != k)
			swap(tmpMatrix[k], tmpMatrix[maxIdx]);

		if (tmpMatrix[k][k] == 0)
			res = Empty;

		for (int i = k + 1; i < rows; i++)
		{
			double coefficient = tmpMatrix[i][k] / tmpMatrix[k][k];
			for (int j = k; j < rows + 1; ++j)
				tmpMatrix[i][j] -= coefficient * tmpMatrix[k][j];
		}
	}

	for (int i = rows - 1; i >= 0; i--) {

		double sum = 0;

		for (int j = i + 1; j < rows; j++)
			sum += tmpMatrix[i][j] * solution[j];
		solution[i] = (tmpMatrix[i][rows] - sum) / tmpMatrix[i][i];
	}

	for (int i = 0; i < rows; i++)
		delete[] tmpMatrix[i];
	delete[] tmpMatrix;

	return res;
}

double calculateDeterminant(double** extendedMatrix, int rows, int cols)
{
	double det = 1.0;

	double** tmpMatrix = new double* [rows];
	for (int i = 0; i < rows; i++)
	{
		tmpMatrix[i] = new double[cols];
		memcpy(tmpMatrix[i], extendedMatrix[i], cols * sizeof(double));
	}

	for (int k = 0; k < rows - 1; ++k)
	{
		double maxValue = std::abs(tmpMatrix[k][k]);
		int maxRow = k;
		int maxCol = k;

		for (unsigned i = k; i < rows; i++)
		{
			for (unsigned j = k; j < rows; j++)
			{
				if (abs(tmpMatrix[i][j]) > maxValue)
				{
					maxValue = abs(tmpMatrix[i][j]);
					maxRow = i;
					maxCol = j;
				}
			}
		}

		if (maxRow != k)
		{
			for (int j = k; j < rows + 1; j++)
				swap(tmpMatrix[k][j], tmpMatrix[maxRow][j]);
			det *= -1.0;
		}

		if (maxCol != k)
		{
			for (int i = k; i < rows; ++i)
				swap(tmpMatrix[i][k], tmpMatrix[i][maxCol]);
			det *= -1.0;
		}

		if (tmpMatrix[k][k] == 0)
		{
			for (int i = 0; i < rows; i++)
				delete[] tmpMatrix[i];
			delete[] tmpMatrix;
			return 0.0;
		}

		for (int i = k + 1; i < rows; ++i)
		{
			double coefficent = tmpMatrix[i][k] / tmpMatrix[k][k];
			for (int j = k; j < rows + 1; ++j)
				tmpMatrix[i][j] -= coefficent * tmpMatrix[k][j];
		}
	}

	for (int i = 0; i < rows; i++)
		det *= tmpMatrix[i][i];

	for (int i = 0; i < rows; i++)
		delete[] tmpMatrix[i];
	delete[] tmpMatrix;

	return det;
}