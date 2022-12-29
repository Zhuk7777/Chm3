#include <iostream>
#include<ctime>
#include"Matrix.h"

double* generateVector(int size, double max, double min);
double scalarMultiplication(double* vec1, double* vec2, int size);
double normOfVec(double* vec, int size);
void dividingVector(double* vec, double value, int size);
double** createHouseholderMatrix(double* w, int size);
double** createMatrixForTest(double** H, double** diag, int size);
double** createDiag(int size);
void printMatrix(double** matr, int size);

int main()
{
	setlocale(0, "rus");
	srand(time(NULL));

	double** diag = createDiag(5);
	double* w = generateVector(5, 13, 5);
	dividingVector(w, normOfVec(w, 5), 5);
	double** H = createHouseholderMatrix(w, 5);
	double** A = createMatrixForTest(H, diag, 5);

	printMatrix(A, 5);
	printMatrix(H, 5);
	printMatrix(diag, 5);

	SymmetricMatrix matr(A, 5);
	int countSpins;
	double* eigenvalues = new double[5];
	double r;

	std::cout << "Введите максимальное допустимое число вращений\n";
	int M;
	std::cin >> M;

	std::cout << "Введите максимальное по модулю значение внедиагональных элементов преобразованной матрицы\n";
	double e;
	std::cin >> e;

	double** T = matr.eigenvectors(e, M, countSpins, eigenvalues, r);

	printMatrix(T, 5);
	for (int i = 0; i < 5; i++)
		std::cout << eigenvalues[i] << " ";
	std::cout << "\n" << r << "\n" << countSpins;


	return 0;




}

double* generateVector(int size, double max, double min) {
	double* x = new double[size];
	for (int j = 0; j < size; j++) {
		x[j] = (rand() % (int)max + min) / (rand() % (int)max + min);
		if (x[j] < min)
			x[j] = x[j] + abs(rand() % (int)max / 3 + min);
	}
	return x;
}


double scalarMultiplication(double* vec1, double* vec2, int size)
{
	double result = 0;
	for (int i = 0; i < size; i++)
		result += vec1[i] * vec2[i];
	return result;
}

double normOfVec(double* vec,int size)
{
	return sqrt(scalarMultiplication(vec, vec, size));

}

void dividingVector(double* vec, double value, int size)
{
	for (int i = 0; i < size; i++)
		vec[i] /= value;
}

double** createHouseholderMatrix(double* w, int size)
{
	double** result = new double* [size];
	for (int i = 0; i < size; i++)
		result[i] = new double[size];

	double** ww = new double* [size];
	for (int i = 0; i < size; i++)
		ww[i] = new double[size];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			ww[i][j] = 2 * w[i] * w[j];

	for(int i=0;i<size;i++)
		for (int j = 0; j < size; j++)
		{
			if (i == j)
				result[i][j] = 1 - ww[i][j];
			else
				result[i][j] = -ww[i][j];
		}

	return result;
}

double** createMatrixForTest(double** H, double** diag, int size)
{
	double** result = new double* [size];
	for (int i = 0; i < size; i++)
		result[i] = new double[size];

	for(int i=0;i<size;i++)
		for (int j = 0; j < size; j++)
		{
			result[i][j] = 0;
			for (int k = 0; k < size; k++)
				if (k == j)
					result[i][j] += H[i][k] * diag[k][j];
		}

	double** result2 = new double* [size];
	for (int i = 0; i < size; i++)
		result2[i] = new double[size];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
		{
			result2[i][j] = 0;
			for (int k = 0; k < size; k++)
				result2[i][j] += result[i][k] * H[j][k];
		}

	return result2;

}

double** createDiag(int size)
{
	double diag[5][5] = { 3.2, 0,    0,    0, 0,
						 0 ,   1.76, 0,    0, 0,
						 0,    0,    1.1,  0, 0,
						 0,    0,    0,   -4, 0,
						 0,    0,    0,    0, 7 };
	double** Diag = new double*[5];
	for (int i = 0; i < size; i++)
		Diag[i] = new double[5];

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			Diag[i][j] = diag[i][j];
	return Diag;
}

void printMatrix(double** matr, int size)
{
	for (int i = 0; i < size; i++)
	{
		std::cout << "\n";
		for (int j = 0; j < size; j++)
			std::cout << matr[i][j] << " ";
	}

	std::cout << "\n";
}