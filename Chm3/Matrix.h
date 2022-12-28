#pragma once
#include <iostream>

using namespace std;

class SymmetricMatrix
{
	double** elements;
	int size;

	void allocateCells(int size)
	{
		elements = new double* [size];
		for (int i = 0; i < size; i++)
			elements[i] = new double[size];

		this->size = size;

	}

	void freeCells()
	{
		for (int i = 0; i < size; i++)
			delete[] elements[i];
		delete[] elements;

		size = 0;
	}

	int sign(double value)
	{
		if (value > 0)
			return 1;
		if (value < 0)
			return -1;
		if (value == 0)
			return 0;
	}


	void calculationCS(int i, int j, double& c, double& s)
	{
		double p = 2 * elements[i][j];
		double q = elements[i][i] - elements[j][j];
		double d = sqrt(p * p + q * q);

		if (q == 0)
		{
			c = sqrt(2) / 2;
			s = sqrt(2) / 2;
		}
		else
		{
			double r = abs(q) / (2 * d);
			c = sqrt(0.5 + r);
			s = sqrt(0.5 - r * sign(p * q));
		}
	}

	double** initializationT0()
	{
		double** T = new double* [size];
		for (int i = 0; i < size; i++)
			T[i] = new double[size];

		for(int i=0;i<size;i++)
			for (int j = 0; j < size; j++)
			{
				if (i == j)
					T[i][j] = 1;
				else
					T[i][j] = 0;
			}

		return T;
	}

	double** initializationA()
	{
		double** A = new double* [size];
		for (int i = 0; i < size; i++)
			A[i] = new double[size];

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				A[i][j] = elements[i][j];

		return A;

	}

	double** multiplicationRightTij(double** matr, double c, double s, int i, int j)
	{
		double** result = new double* [size];
		for (int i = 0; i < size; i++)
			result[i] = new double[size];

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				result[i][j] = matr[i][j];

		for (int l = 0; l < size; l++)
		{
			result[l][i] = c * matr[l][i] + s * matr[l][j];
			result[l][j] = -s * matr[l][i] + c * matr[l][j];
		}
	}

	double** multiplicationLeftTijTransposed(double** matr, double c, double s, int i, int j)
	{
		double** result = new double* [size];
		for (int i = 0; i < size; i++)
			result[i] = new double[size];

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				result[i][j] = matr[i][j];

		for (int m = 0; m < size; m++)
		{
			result[i][m] = c * matr[i][m] + s * matr[j][m];
			result[j][m] = -s * matr[i][m] + c * matr[j][m];
		}

		return result;
	}
public:

	SymmetricMatrix()
	{
		size = 0;
		elements = nullptr;
	}

	SymmetricMatrix(const const SymmetricMatrix& M)
	{
		allocateCells(M.size);
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				elements[i][j] = M.elements[i][j];
	}

	SymmetricMatrix(int size,double* list)
	{
		allocateCells(size);
		int k = 0;

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				elements[i][j] = list[k];
				k++;
			}
	}

	~SymmetricMatrix()
	{
		freeCells();
	}

	void setMatrix(int size, double* list)
	{
		allocateCells(size);
		int k = 0;

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				elements[i][j] = list[k];
				k++;
			}
	}

	void fillMatrix()
	{
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				std::cin >> elements[i][j];
	}

	double absMax(int& _i, int& _j)
	{
		double max = -1;
		for (int i = 0; i < size; i++)
			for (int j = 0; j < i; j++)
				if (abs(elements[i][j]) > max)
				{
					max = abs(elements[i][j]);
					_i = i;
					_j = j;
				}

		return max;
	}

	double** eigenvectors(double e, int M, int& K, double* a, double& r)
	{
		double** T = initializationT0();
		double** A = initializationA();
		int iterator = 0;
		int i, j;
		double E = absMax(i, j);

		while (iterator<M && e<E)
		{
		}
	}
};