#pragma once
#include <iostream>

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
			s = sqrt(0.5 - r)*sign(p * q);
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

		return result;
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

	void creatingVectorOfEigenvalues(double** A, double* a)
	{
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				if (i == j)
					a[i] = A[i][j];


	}

	double measureOfAccuracy(double**A,double**T)
	{
		double** result = new double*[size];
		for (int i = 0; i < size; i++)
			result[i] = new double[size];

		double** leftPart = new double* [size];
		for (int i = 0; i < size; i++)
			leftPart[i] = new double[size];

		double** rightPart = new double* [size];
		for (int i = 0; i < size; i++)
			rightPart[i] = new double[size];

		double max = -1;

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
			{
				leftPart[i][j] = T[i][j] * A[j][j];
				if (j <= i)
				{
					rightPart[i][j] = 0;
					for (int k = 0; k < size; k++)
						rightPart[i][j] += elements[i][k] * T[k][j];

					rightPart[j][i] = rightPart[i][j];
				}

				result[i][j] = leftPart[i][j] - rightPart[i][j];
				if (abs(result[i][j]) > max)
					max = abs(result[i][j]);
			}

		return max;
	
	}

	void print(double** matr)
	{
		for (int i = 0; i < size; i++)
		{
			std::cout << "\n";
			for (int j = 0; j < size; j++)
				std::cout << matr[i][j] << " ";
		}
		std::cout << "\n";
	}

	double absMax(double** matr, int& _i, int& _j)
	{
		double max = -1;
		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				if (i!=j&&abs(matr[i][j]) > max)
				{
					max = abs(matr[i][j]);
					_i = i;
					_j = j;
				}

		return max;
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

	SymmetricMatrix(double** matr, int size)
	{
		allocateCells(size);

		for (int i = 0; i < size; i++)
			for (int j = 0; j < size; j++)
				elements[i][j] = matr[i][j];
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

	void printMatrix()
	{
		for (int i = 0; i < size; i++)
		{
			std::cout << "\n";
			for (int j = 0; j < size; j++)
				std::cout << elements[i][j] << " ";
		}
		std::cout << "\n";

	}

	double** eigenvectors(double e, int M, int& K, double* a, double& r)
	{
		double** T = initializationT0();
		double** A = initializationA();
		int iterator = 0;
		int i, j;
		double c, s;
		double E = absMax(A, i, j);

		while (iterator<M && e<E)
		{
			E = absMax(A, i, j);
			calculationCS(i, j, c, s);
			A = multiplicationLeftTijTransposed(A, c, s, i, j);
			//print(A);
			A = multiplicationRightTij(A, c, s, i, j);
			T = multiplicationRightTij(T, c, s, i, j);
			iterator++;
			//print(A);

		}

		K = iterator;
		creatingVectorOfEigenvalues(A, a);
		r = measureOfAccuracy(A, T);

		return T;
	}
};