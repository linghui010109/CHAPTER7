#include <iostream>
#include<cmath>
#include<vector>
#include <string>
#include"Eigen/Dense"
#include"Function.h"
using namespace std;
using namespace Eigen;

class SolveMatrix
{
public:
	MatrixXd solve_matrix(MatrixXd A, MatrixXd b) {
	MatrixXd x;
	x = A.fullPivLu().solve(b);
	return x;
	}
};

class Error
{
public:
	Error(MatrixXd A, double* B, int norm)
	{
		int m = A.cols();
		double e = 0;
		if (norm == 1)
		{
			for (int i = 0; i < m; i++)
			{
				e += fabs( A(i,0) - B[i]);
			}
		}
		else if (norm == 2)
		{
			for (int i = 0; i < m; i++)
			{
				e += (A(i, 0) - B[i]) * (A(i, 0) - B[i]);
			}
			e = sqrt(e);
		}
		else
		{
			for (int i = 0; i < m; i++)
			{
				double d = fabs(A(i, 0) - B[i]);
				if (e < d)
				{
					e = d;
				}
			}
		}

		
	}
};

class PDE
{
public:
	double** set_left_matrix(int n);
	double** Set_right_matrix(int n, Function& fun1);
	PDE(int n, Function& func, Function& funcdiff, const string& type_name);
	MatrixXd transfer_A(double** A,int n);
	MatrixXd transfer_F(double** F,int n);
protected:
	int n;
	double** A;
	double** F;
	int s;
	int m;
	double h;
};

MatrixXd PDE::transfer_A(double** A,int n)
{
	MatrixXd A2D((n-1)*(n-1), (n-1)*(n-1));

	for (int i = 0; i < (n - 1) * (n - 1); i++)
	{
		for (int j = 0; j < (n - 1) * (n - 1); j++)
		{
			A2D(i, j) = A[i+1][j+1];
			cout << A2D(i, j) << " ";
		}
		cout << endl;
	}
	return A2D;
}

MatrixXd PDE::transfer_F(double** F,int n)
{
	//int Alen = sizeof(F) / sizeof(double);
	//int Alen2 = sizeof(F[0]) / sizeof(double);//第二维数组的长度
	//int Alen3 = Alen / Alen2;//第一维数组的长度
	//cout <<"Alen:"<< Alen << endl;
	//cout << "Alen2:" << Alen2 << endl;
	//cout << "Alen3:" << Alen3 << endl;

	MatrixXd F2D((n-1)*(n-1), 1);
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			F2D((i-1) + (j-1)* (n-1) , 0) = F[i][j];
		}
	}
	cout << "F:" << endl;
	for (int i = 0; i < (n - 1) * (n - 1); i++)
	{
		cout << F2D(i, 0) << endl;
	}
	return F2D;

}

//矩阵A的内部
double** PDE::set_left_matrix(int n)
{
	int m = n - 1;
	int s = (n + 1) * (n + 1);
	double** A;
	A = new double* [s];
	double h = 1.0 / n;
	for (int k = 0; k < s; k++)
	{
		A[k] = new double[s];
	}
	for (int i = 1; i <= (n - 1) * (n - 1); i++)
	{
		for (int j = 1; j <= i; j++)
		{
			if (i == j)
			{
				A[i][j] = 4.0;
			}
			else if (j == i - 3)
			{
				A[i][j] = -1.0;
			}
			else if (i == j + 1)
			{
				if ((j) % (n-1) == 0)
				{
					A[i][j] = 0;
				}
				else { A[i][j] = -1.0; }
			}
			else { A[i][j] = 0; }
			A[j][i] = A[i][j];
		}
	}
	return A;
	//A[0][0] = 0;
	/*for (int i = 0; i < s; i++)
	{
		delete[] A[i];
	}
	delete[] A;*/
	/*for (int i = 1; i <= (n - 1) * (n - 1); i++)
	{
		for (int j = 1; j <= (n-1)*(n-1); j++)
		{

			cout << A[i][j]<<" ";
		}
		cout << endl;
	}*/
	//cout << A[0][0];
}

//矩阵F的初始化
double** PDE::Set_right_matrix(int n, Function& f1)//偏微分方程
{
	//SolveMatrix s1;
	int s = (n + 1) * (n + 1);
	double** F;
	F = new double* [s];
	double h = 1.0 / n;
	for (int k = 0; k < s; k++)
	{
		F[k] = new double[s];
	}
	for (int i = 1; i <= (n - 1); i++)
	{
		for (int j = 1; j <= (n - 1); j++)
		{
			F[i][j] =h*h*f1(i * h, j * h);
			//cout << h*h*F[i][j]<<" ";
		}
		//cout << endl;
	}
	return F;
}

PDE::PDE(int n, Function& func, Function& funcdiff, const string& type_name)//func为原方程
{
	SolveMatrix s1;
	double** A = set_left_matrix(n);
	double** F = Set_right_matrix(n, funcdiff);
	double h = 1.0 / n;
	int m = n - 1;
	/*for (int i = 1; i <= (n - 1); i++)
	{
		for (int j = 1; j <= (n - 1); j++)
		{
			cout << F[i][j] << " ";
		}
		cout << endl;
	}*/
	if (type_name == "Dirichlet")
	{
		//cout << F[1][1];
		//cout << func(0.0, h);

		//四个角落
		F[1][1] = F[1][1] + func(0.0, h) + func(h, 0.0);
		F[m][m] = F[m][m] + func(1.0, m * h) + func(m * h, 1.0);
		F[1][m] = F[1][m] + func(h, 1.0) + func(0.0, m*h);
		F[m][1] = F[m][1] + func(1.0,h) + func(m * h,0.0);
		//cout << F[1][m] << endl;
		//最下边
		for (int i = 2; i < m; i++)
		{
			F[i][1] = F[i][1] + func(i * h, 0.0);
		}

		//最右边
		for (int i = 2; i < m; i++)
		{
			F[m][i] = F[m][i] + func(1.0, i * h);
		}

		//最左边
		for (int i = 2; i < m; i++)
		{
			F[1][i] = F[1][i] + func(0.0, i * h);
		}

		//最上边
		for (int i = 2; i < m; i++)
		{
			F[i][m] = F[i][m] + func(i * h, 1.0);
		}
		/*cout << "new" << endl;
		for (int i = 1; i <= (n - 1); i++)
		{
			for (int j = 1; j <= (n - 1); j++)
			{
				cout << F[i][j] << " ";
			}
			cout << endl;
		}*/
		//cout << "A:" << endl;
		MatrixXd A2D = transfer_A(A, n);
		MatrixXd F2D = transfer_F(F,n);
		MatrixXd u=s1.solve_matrix(A2D, F2D);
		cout << "u:" << endl;
		cout << u << endl;
		cal_true_solution c1(4, func);
		//Error err(u,);
	}
}

class cal_true_solution
{
public:
	cal_true_solution(int n, Function& fun)
	{
		double h = 1.0 / n;
		double* real_u = new double[(n-1)*(n-1)];
		cout << "real u:" << endl;
		for (int i = 0; i < n-1; i++)
		{
			for (int j = 0; j < n-1; j++)
			{
				real_u[ i + j  * (n - 1)] = fun((i+1) * h, (j+1) * h);
				
				cout << real_u[i + j * (n - 1)] << endl;
			}
			
		}

	}
};




