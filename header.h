#include <iostream>
#include<cmath>
#include<vector>
#include"Eigen/Dense"
#include"Function.h"
using namespace std;
using namespace Eigen;


class PDE
{
public:
	double** set_left_matrix(int n);
	double** Set_right_matrix(int n, Function& fun1);
	PDE(int n, Function& func, Function& funcdiff, const string& type_name);
protected:
	int n;
	double** A;
	double** F;
	int s;
	int m;
	double h;
};

//矩阵A的内部
double** PDE::set_left_matrix(int n)
{
	int m = n - 1;
	int s = (n+1)*(n+1);
	double** A;
	A = new double* [s];
	double h = 1.0 / n;
	for (int k = 0; k < s; k++)
	{
		A[k] = new double [s];
	}
	for (int i = 1; i <= (n-1)*(n-1) ; i++)
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
				if ((j ) % n == 0)
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
double** PDE::Set_right_matrix(int n,Function& f1)//偏微分方程
{
	int s = (n + 1) * (n + 1);
	double** F;
	F = new double* [s];
	double h = 1.0 / n;
	for (int k = 0; k < s; k++)
	{
		F[k] = new double[s];
	}
	for (int i = 1; i <= (n-1); i++)
	{
		for (int j = 1; j <= (n-1); j++)
		{
			F[i][j] = f1(i * h, j * h);
			//cout << h*h*F[i][j]<<" ";
		}
		//cout << endl;
	}
	return F;
}

PDE::PDE(int n, Function& func, Function& funcdiff,const string& type_name)//func为原方程
{
	double** A=set_left_matrix(n);
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
	if (type_name == "Direchlet")
	{
		cout << F[1][1];
		cout << func(0.0, h);

		//四个角落
		F[1][1] = F[1][1] + func(0.0, h) + func(h, 0.0);

		F[m][m] = F[m][m] + func(1.0, m*h)+func(m*h,1.0);
		F[1][m] = F[1][m] + func(1.0, h) + func(m*h, 0.0);
		F[m][1] = F[m][1] + func(h, 1.0) + func(0.0, m * h);
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
		for (int i=2;i<m;i++)
		{
			F[i][m] = F[i][m] + func(i * h, 1.0);
		}

		for (int i = 1; i <= (n - 1); i++)
		{
			for (int j = 1; j <= (n - 1); j++)
			{
				cout <<F[i][j]<<" ";
			}
			cout << endl;
		}

	}
}

//class pde 
//{
//public:
//	void getvalue(int n,function& func,const string& type_name)
//	{
//		if (type_name == "direchlet")
//		{
//			double h = 1.0 / n;
//			cout << func(0.0, h);
//			cout << f[1][1];
//			//四个角落
//			f[1][1] = f[1][1] + func(0.0, h) + func(h, 0.0);
//			f[m][m] = f[m][m] + func(1.0, m*h)+func(m*h,1.0);
//			f[1][m] = f[1][m] + func(1.0, h) + func(m*h, 0.0);
//			f[m][1] = f[m][1] + func(h, 1.0) + func(0.0, m * h);
//
//			//最下边
//			for (int i = 2; i < m; i++)
//			{
//				f[i][1] = f[i][1] + func(i * h, 0.0);
//			}
//
//			//最右边
//			for (int i = 2; i < m; i++)
//			{
//				f[m][i] = f[m][i] + func(1.0, i * h);
//			}
//
//			//最左边
//			for (int i = 2; i < m; i++)
//			{
//				f[1][i] = f[1][i] + func(0.0, i * h);
//			}
//			
//			//最上边
//			for (int i=2;i<m;i++)
//			{
//				f[i][m] = f[i][m] + func(i * h, 1.0);
//			}
//		}
//	}
//};

class Solver
{
public:
	virtual void Dirichlet_Solver(int h,int n,MatrixXf& left,MatrixXf& right);
	//void Nuemann_Solver();
	//void Mix_Solver();

};
void Solver::Dirichlet_Solver(int h, int n, MatrixXf& left, MatrixXf& right)
{

}
