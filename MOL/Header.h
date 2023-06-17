#include<iostream>
#include<cmath>
#include"Eigen/Dense"
//s#include"Function.h"

using namespace std;
using namespace Eigen;
class A
{
public:
	MatrixXd create_A(int n, double epsilon, double h,double a);
}; 

MatrixXd A::create_A(int n, double epsilon, double h,double a)
{
	MatrixXd A_l(n - 1, n - 1);
	A_l = MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = 0; j < n - 1; j++)
		{
			if (i = j + 1)
			{
				A_l(i, j) = a / (2 * h);
			}
			A_l(j, i) = -A_l(i, j);
			A_l(0, n - 1) = a / (2 * h);
			A_l(n - 1, 0) = -a / (2 * h);
		}
	}
	MatrixXd A_r(n-1, n-1);
	A_r = MatrixXd::Zero(n, n);
	A_r = (-2*epsilon / pow(h, 2)) * MatrixXd::Identity(n, n);
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = 0; j < n - 1; j++)
		{
			if (i = j + 1)
			{
				A_r(i, j) = epsilon / pow(h, 2);
			}
			A_r(n - 1, 0) = epsilon / (h * h);
			A_r(0, n - 1) = epsilon / (h * h);
			A_r(j, i) = A_r(i, j);
		}
	}
	cout << "A_r" << endl;
	cout << A_r << endl;
	cout << "A_l" << endl;
	cout << A_l << endl;
	return A_r + A_l;
}

class Lax_Friedrichs
{
public:
	MatrixXd solve_lax_Friedrichs(int n, double a, double h, double k, int T, MatrixXd u_0);
};

MatrixXd Lax_Friedrichs::solve_lax_Friedrichs(int n, double a, double h, double k, int T, MatrixXd u_0)
{
	A a1;
	double epsilon = h * h / (2 * k);
	//cout << epsilon << endl;
	MatrixXd A_lax = a1.create_A(n, epsilon, h, a);
	MatrixXd R = MatrixXd::Identity(n, n) + A_lax * k;
	MatrixXd u = R * u_0;
	return u;
}

class Lax_Wendroff
{
public:
	MatrixXd solve_lax_Wendroff(int n, double a, double h, double k, int T, MatrixXd u_0);
};

MatrixXd Lax_Wendroff::solve_lax_Wendroff(int n, double a, double h, double k, int T, MatrixXd u_0)
{
	A a1;
	double epsilon = k*a*a/2;
	//cout << epsilon << endl;
	MatrixXd A_lax = a1.create_A(n, epsilon, h, a);
	MatrixXd R = MatrixXd::Identity(n, n) + A_lax * k;
	MatrixXd u = R * u_0;
	return u;
}

class upwind
{
public:
	MatrixXd solve_upwind(int n, double a, double h, double k, int T, MatrixXd u_0);
};

MatrixXd upwind::solve_upwind(int n, double a, double h, double k, int T, MatrixXd u_0)
{
	A a1;
	double epsilon = h * a / 2;
	//cout << epsilon << endl;
	MatrixXd A_lax = a1.create_A(n, epsilon, h, a);
	MatrixXd R = MatrixXd::Identity(n, n) + A_lax * k;
	MatrixXd u = R * u_0;
	return u;
}
