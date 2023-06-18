#include<iostream>
#include<cmath>
#include"Eigen/Dense"
#include"Function.h"

using namespace std;
using namespace Eigen;
class A
{
public:
	MatrixXd create_A(int n, double v, double h);
};

MatrixXd A::create_A(int n, double v, double h)
{
	MatrixXd A(n, n);
	A = (-2 * v / pow(h, 2)) * MatrixXd::Identity(n, n);
	for (int i = 0; i < n-1; i++)
	{
		for (int j = 0; j < n-1; j++)
		{
			if (i = j + 1)
			{
				A(i, j) = v / pow(h, 2);
			}
			
			A(j, i) = A(i, j);
		}
	}
	//cout << A;
	return A;
}

class Crank_Nicolson
{
public:
	MatrixXd solve_crank_nicolson(int n, double v, double h,double r,MatrixXd u_0,int t);
};

MatrixXd Crank_Nicolson::solve_crank_nicolson(int n, double v, double h,double r,MatrixXd u_0,int t)
{
	A a1;
	double k = r * pow(h, 2) / v;
	MatrixXd a = a1.create_A(n-1, v, h);
	//cout << "a" << endl;
	//cout << a << endl;
	MatrixXd l = MatrixXd::Identity(n-1, n-1) - a * k / 2.0;
	//cout << "l" << endl;
	//cout << l << endl;
	MatrixXd R = MatrixXd::Identity(n-1, n-1) + a * k / 2.0;
	//cout << "R" << endl;
	//cout << R << endl;
	MatrixXd u = l.inverse() * R * u_0;
	//for (int i = 1; i < t; i++)
	if(t>1)
	{
		Crank_Nicolson t1;
		MatrixXd uu=t1.solve_crank_nicolson(n, v, h, r, u, t-1);
		//cout << "uu" << i << uu << endl;
		u = uu;
	}
	//cout << "u" << endl;
	//cout << u << endl;

	return u;
}

class FTCS
{
public:
	MatrixXd solve_FTCS(int n, double v, double h, double r, MatrixXd u_0,int t);
};

MatrixXd FTCS::solve_FTCS(int n, double v, double h, double r, MatrixXd u_0,int t)
{
	A a1;
	double k = r * pow(h, 2) / v;
	MatrixXd a = a1.create_A(n - 1, v, h);
	MatrixXd R = MatrixXd::Identity(n - 1, n - 1) + a * k ;
	MatrixXd u =  R * u_0;
	if (t > 1)
	{
		FTCS f1;
		MatrixXd uu = f1.solve_FTCS(n, v, h, r, u, t - 1);
		//cout << "uu" << i << uu << endl;
		u = uu;
	}
	//cout << "u" << endl;
	//cout << u << endl;
	return u;
}

class BTCS
{
public:
	MatrixXd solve_BTCS(int n, double v, double h, double r, MatrixXd u_0, int t);
};

MatrixXd BTCS::solve_BTCS(int n, double v, double h, double r, MatrixXd u_0, int t)
{
	A a1;
	double k = r * pow(h, 2) / v;
	MatrixXd a = a1.create_A(n - 1, v, h);
	MatrixXd l = MatrixXd::Identity(n - 1, n - 1) - a * k ;
	MatrixXd u = l.inverse() * u_0;
	if (t > 1)
	{
		BTCS f1;
		MatrixXd uu = f1.solve_BTCS(n, v, h, r, u, t - 1);
		//cout << "uu" << i << uu << endl;
		u = uu;
	}
	//cout << "u" << endl;
	//cout << u << endl;
	return u;
}
