#include <iostream>
#include <cmath>
#include <vector>
#include"Eigen/Dense"
#include "Adam_bashforth.h"
#define Epsilon 0.000000001

using namespace std;
using namespace Eigen;

class Newton_Iter
{
public:
	MatrixXd ff(vector<double> u, vector<double> u_newton, double k);//原函数p=2
	MatrixXd Jacobian(vector<double> u, double k, int p);//jacobi矩阵
	vector<double> Newton(vector<double> newton_initial_value, vector<double> u, vector<double> u2, vector<double> u3, vector<double> u4, int n, double k, int p);//n为迭代次数
	MatrixXd f3(vector<double> u1, vector<double> u2, vector<double> u_newton, double k);//原函数p=3
	MatrixXd f4(vector<double> u1, vector<double> u2, vector<double> u3, vector<double> u_newton, double k);//原函数p=4
	MatrixXd f5(vector<double> u1, vector<double> u2, vector<double> u3, vector<double>, vector<double> u_newton, double k);//原函数p=5
};

MatrixXd Newton_Iter::ff(vector<double> u, vector<double> u_newton, double k)
{
	double mu = 0.01227747;
	vector<double> v;
	MatrixXd vv(6, 1);
	vector<double> value = f(u, 0);
	vector<double> newton_value = f(u_newton, 0);

	v = u_newton - u - 0.5 * k * value - 0.5 * k * newton_value;
	for (int i = 0; i < 6; i++)
	{
		vv(i, 0) = v[i];
	}
	//cout << "vv:" << vv << endl;
	return vv;
}

MatrixXd Newton_Iter::f3(vector<double> u1, vector<double> u2, vector<double> u_newton, double k)
{
	double mu = 0.01227747;
	vector<double> v;
	MatrixXd vv(6, 1);
	vector<double> value1 = f(u1, 0);
	vector<double> value2 = f(u2, 0);
	vector<double> newton_value = f(u_newton, 0);

	v = u_newton - value2 + (1.0 / 12.0) * k * value1 - (8.0 / 12.0) * k * value2 - (5.0 / 12.0) * k * newton_value;
	for (int i = 0; i < 6; i++)
	{
		vv(i, 0) = v[i];
	}
	//cout << "vv:" << vv << endl;
	return vv;
}

MatrixXd Newton_Iter::f4(vector<double> u1, vector<double> u2, vector<double> u3, vector<double> u_newton, double k)
{
	double mu = 0.01227747;
	vector<double> v;
	MatrixXd vv(6, 1);
	vector<double> value1 = f(u1, 0);
	vector<double> value2 = f(u2, 0);
	vector<double> value3 = f(u3, 0);
	vector<double> newton_value = f(u_newton, 0);

	v = u_newton - value3 - (1.0 / 24.0) * k * value1 + (5.0 / 24.0) * k * value2 - (19.0 / 24.0) * k * value3 - (9.0 / 24.0) * k * newton_value;
	for (int i = 0; i < 6; i++)
	{
		vv(i, 0) = v[i];
	}
	//cout << "vv:" << vv << endl;
	return vv;
}

MatrixXd Newton_Iter::f5(vector<double> u1, vector<double> u2, vector<double> u3, vector<double> u4, vector<double> u_newton, double k)
{
	double mu = 0.01227747;
	vector<double> v;
	MatrixXd vv(6, 1);
	vector<double> value1 = f(u1, 0);
	vector<double> value2 = f(u2, 0);
	vector<double> value3 = f(u3, 0);
	vector<double> value4 = f(u4, 0);
	vector<double> newton_value = f(u_newton, 0);

	v = u_newton - value4 + (19.0 / 720.0) * k * value1 - (106.0 / 720.0) * k * value2 + (264.0 / 720.0) * k * value3 - (646.0 / 720.0) * k * value4 - (251.0 / 720.0) * k * newton_value;
	for (int i = 0; i < 6; i++)
	{
		vv(i, 0) = v[i];
	}
	//cout << "vv:" << vv << endl;
	return vv;
}


MatrixXd Newton_Iter::Jacobian(vector<double> u, double k, int p)
{

	double mu = 0.01227747;
	MatrixXd j(6, 6);
	if (p == 2)
	{
		j(0, 0) = 1;
		j(1, 0) = j(2, 0) = 0;
		j(3, 0) = -0.5 * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) + (3 * mu * (2 * mu + 2 * u[0] - 2) * (mu + u[0] - 1)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * (2 * mu + 2 * u[0]) * (mu + u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) + 1);
		j(4, 0) = -0.5 * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(5, 0) = -0.5 * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(0, 1) = j(2, 1) = 0;
		j(1, 1) = j(2, 2) = j(3, 3) = j(4, 4) = j(5, 5) = 1;
		j(3, 1) = -0.5 * k * ((3 * mu * u[1] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 1) = -0.5 * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[1], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[1], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + 1);
		j(5, 1) = -0.5 * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(0, 2) = j(1, 2) = 0;
		j(3, 2) = -0.5 * k * ((3 * mu * u[2] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[2] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 2) = -0.5 * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(5, 2) = -0.5 * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[2], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[2], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(1, 3) = j(2, 3) = j(5, 3) = 0;
		j(4, 3) = 2.0 * 0.5 * k;
		j(3, 4) = -2.0 * 0.5 * k;
		j(0, 3) = j(1, 4) = j(2, 5) = -0.5 * k;
		j(0, 4) = j(2, 4) = j(5, 4) = 0;
		j(0, 5) = j(1, 5) = j(3, 5) = j(4, 5) = 0;
		//cout << "j:" << endl;
		//cout << j << endl;
		MatrixXd inverse_j = j.inverse();
		//cout << "j_inverse:" << endl;
		//cout << inverse_j << endl;
		return inverse_j;
		//return j;
	}
	else if (p == 3)
	{
		j(0, 0) = 1;
		j(1, 0) = j(2, 0) = 0;
		j(3, 0) = -(5.0 / 12.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) + (3 * mu * (2 * mu + 2 * u[0] - 2) * (mu + u[0] - 1)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * (2 * mu + 2 * u[0]) * (mu + u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) + 1);
		j(4, 0) = -(5.0 / 12.0) * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(5, 0) = -(5.0 / 12.0) * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(0, 1) = j(2, 1) = 0;
		j(1, 1) = j(2, 2) = j(3, 3) = j(4, 4) = j(5, 5) = 1;
		j(3, 1) = -(5.0 / 12.0) * k * ((3 * mu * u[1] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 1) = -(5.0 / 12.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[1], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[1], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + 1);
		j(5, 1) = -(5.0 / 12.0) * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(0, 2) = j(1, 2) = 0;
		j(3, 2) = -(5.0 / 12.0) * k * ((3 * mu * u[2] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[2] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 2) = -(5.0 / 12.0) * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(5, 2) = -(5.0 / 12.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[2], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[2], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(1, 3) = j(2, 3) = j(5, 3) = 0;
		j(4, 3) = 2.0 * (5.0 / 12.0) * k;
		j(3, 4) = -2.0 * (5.0 / 12.0) * k;
		j(0, 3) = j(1, 4) = j(2, 5) = -(5.0 / 12.0) * k;
		j(0, 4) = j(2, 4) = j(5, 4) = 0;
		j(0, 5) = j(1, 5) = j(3, 5) = j(4, 5) = 0;
		MatrixXd inverse_j = j.inverse();
		return inverse_j;
	}

	else if (p == 4)
	{
		j(0, 0) = 1;
		j(1, 0) = j(2, 0) = 0;
		j(3, 0) = -(9.0 / 24.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) + (3 * mu * (2 * mu + 2 * u[0] - 2) * (mu + u[0] - 1)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * (2 * mu + 2 * u[0]) * (mu + u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) + 1);
		j(4, 0) = -(9.0 / 24.0) * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(5, 0) = -(9.0 / 24.0) * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(0, 1) = j(2, 1) = 0;
		j(1, 1) = j(2, 2) = j(3, 3) = j(4, 4) = j(5, 5) = 1;
		j(3, 1) = -(9.0 / 24.0) * k * ((3 * mu * u[1] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 1) = -(9.0 / 24.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[1], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[1], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + 1);
		j(5, 1) = -(9.0 / 24.0) * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(0, 2) = j(1, 2) = 0;
		j(3, 2) = -(9.0 / 24.0) * k * ((3 * mu * u[2] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[2] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 2) = -(9.0 / 24.0) * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(5, 2) = -(9.0 / 24.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[2], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[2], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(1, 3) = j(2, 3) = j(5, 3) = 0;
		j(4, 3) = 2.0 * (9.0 / 24.0) * k;
		j(3, 4) = -2.0 * (9.0 / 24.0) * k;
		j(0, 3) = j(1, 4) = j(2, 5) = -(9.0 / 24.0) * k;
		j(0, 4) = j(2, 4) = j(5, 4) = 0;
		j(0, 5) = j(1, 5) = j(3, 5) = j(4, 5) = 0;
		MatrixXd inverse_j = j.inverse();
		return inverse_j;
	}
	else if (p == 5)
	{
		j(0, 0) = 1;
		j(1, 0) = j(2, 0) = 0;
		j(3, 0) = -(251.0 / 720.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) + (3 * mu * (2 * mu + 2 * u[0] - 2) * (mu + u[0] - 1)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * (2 * mu + 2 * u[0]) * (mu + u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) + 1);
		j(4, 0) = -(251.0 / 720.0) * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(5, 0) = -(251.0 / 720.0) * k * ((3 * mu * u[1] * (2 * mu + 2 * u[0] - 2)) / (2 * pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0))) - (3 * u[1] * (2 * mu + 2 * u[0]) * (mu - 1)) / (2 * pow((pow((mu + u[0]), 2) + pow(u[2], 2) + pow(u[2], 2)), (5.0 / 2.0))));
		j(0, 1) = j(2, 1) = 0;
		j(1, 1) = j(2, 2) = j(3, 3) = j(4, 4) = j(5, 5) = 1;
		j(3, 1) = -(251.0 / 720.0) * k * ((3 * mu * u[1] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 1) = -(251.0 / 720.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[1], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[1], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + 1);
		j(5, 1) = -(251.0 / 720.0) * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(0, 2) = j(1, 2) = 0;
		j(3, 2) = -(251.0 / 720.0) * k * ((3 * mu * u[2] * (mu + u[0] - 1)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[2] * (mu + u[0]) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(4, 2) = -(251.0 / 720.0) * k * ((3 * mu * u[1] * u[2]) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) - (3 * u[1] * u[2] * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(5, 2) = -(251.0 / 720.0) * k * ((mu - 1) / pow((pow(mu + u[0], 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - mu / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (3.0 / 2.0)) - (3 * pow(u[2], 2) * (mu - 1)) / pow((pow((mu + u[0]), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)) + (3 * mu * pow(u[2], 2)) / pow((pow((mu + u[0] - 1), 2) + pow(u[1], 2) + pow(u[2], 2)), (5.0 / 2.0)));
		j(1, 3) = j(2, 3) = j(5, 3) = 0;
		j(4, 3) = 2.0 * (251.0 / 720.0) * k;
		j(3, 4) = -2.0 * (251.0 / 720.0) * k;
		j(0, 3) = j(1, 4) = j(2, 5) = -(251.0 / 720.0) * k;
		j(0, 4) = j(2, 4) = j(5, 4) = 0;
		j(0, 5) = j(1, 5) = j(3, 5) = j(4, 5) = 0;
		MatrixXd inverse_j = j.inverse();
		return inverse_j;
	}

}

vector<double> Newton_Iter::Newton(vector<double> newton_initial_value, vector<double> u, vector<double> u2, vector<double> u3, vector<double> u4, int n, double k, int p)
{
	int iter = 0;
	vector<double> x(6);

	do
	{
		iter = iter + 1;
		MatrixXd f_vec(6, 1);
		MatrixXd jacobian_mat(6, 6);
		if (p == 2)
		{
			f_vec = ff(u, newton_initial_value, k);
		}
		else if (p == 3)
		{
			f_vec = f3(u, u2, newton_initial_value, k);
		}
		else if (p == 4)
		{
			f_vec = f4(u, u2, u3, newton_initial_value, k);
		}
		else if (p == 5)
		{
			f_vec = f5(u, u2, u3, u4, newton_initial_value, k);
		}
		//cout << f_vec << endl;
		jacobian_mat = Jacobian(newton_initial_value, k, p);
		//cout << "j_inverse:" << endl;
		//cout << jacobian_mat << endl;
		vector<double> u_new(6);
		double sum = 0;
		for (int i = 0; i < 6; i++)
		{
			sum = 0;
			for (int j = 0; j < 6; j++)
			{
				sum = sum + jacobian_mat(i, j) * f_vec(j);
			}
			u_new[i] = newton_initial_value[i] - sum;
		}
		double errornorm = 0;
		for (int i = 0; i < 6; i++)
		{
			errornorm = errornorm + fabs(u_new[i] - newton_initial_value[i]);
		}
		if (errornorm < Epsilon)
		{
			break;
		}
		for (int i = 0; i < 6; i++)
		{
			newton_initial_value[i] = u_new[i];
		}
	} while (iter < n);
	/*for (int i = 0; i < 6; i++)
	{
		cout << u0[i] << endl;
	}*/
	return newton_initial_value;
}

class Adam_Moulton
{
private:
	int n;
	double T, mu;
public:
	vector<double> Adam_Moulton_solver(int n, vector<double>& u, double T, int p);
};

vector<double> Adam_Moulton::Adam_Moulton_solver(int n, vector<double>& u, double T, int p)
{
	classical_RK r2;
	Newton_Iter i1;
	double k = 1.0 * T / n;
	
	if (p == 2)
	{
		vector<double> u2, u3, u4;
		vector<double> u_new;
		vector<double> init;
		
		for (int i = 1; i < n; i++)
		{
			vector<double> init;
			for (int j = 0; j < 6; j++)
			{
				init.push_back(u[j]);
			}
			u_new=i1.Newton(init,u,u2,u3,u4, 5,k,2);
			u = u_new;
			cout << i << ":";
			for (int j = 0; j < 6; j++)
			{
				cout << u[j]<<" ";
			}
			cout << endl;
		}
		return u;
	}
	else if (p == 3)
	{
		vector<double> u1 = r2.classical_RK_solver(n, u, T, 1);
		vector<double> u_new,u3,u4;
		
		for (int i = 2; i < n; i++)
		{	
			vector<double> init;
			for (int j = 0; j < 6; j++)
			{
				init.push_back(u1[j]);
			}
			u_new = i1.Newton(init, u, u1, u3, u4, 5, k, 3);
			u = u1;
			u1 = u_new;

			cout << i << ":";
			for (int j = 0; j < 6; j++)
			{
				cout << u1[j] << " ";
			}
			cout << endl;
		}
		return u1;
	}
	else if (p == 4)
	{
		vector<double> u1 = r2.classical_RK_solver(n, u, T, 1);
		vector<double> u2 = r2.classical_RK_solver(n, u, T, 2);
		vector<double> u_new,u4;

		for (int i = 3; i < n; i++)
		{
			vector<double> init;
			for (int j = 0; j < 6; j++)
			{
				init.push_back(u2[j]);
			}
			u_new = i1.Newton(init, u, u1, u2, u4, 5, k, 3);
			u = u1;
			u1 = u2;
			u2 = u_new;
			cout << i << ":";
			for (int j = 0; j < 6; j++)
			{
				cout << u2[j] << " ";
			}
			cout << endl;
		}
		return u2;
	}
	else if (p == 5)
	{
		vector<double> u1 = r2.classical_RK_solver(n, u, T, 1);
		vector<double> u2 = r2.classical_RK_solver(n, u, T, 2);
		vector<double> u3 = r2.classical_RK_solver(n, u, T, 5);
		vector<double> u_new;
		
		for (int i = 4; i < n; i++)
		{
			vector<double> init;
			for (int j = 0; j < 6; j++)
			{
				init.push_back(u3[j]);
			}
			u_new = i1.Newton(init, u, u1, u2, u3, 5, k,5);
			u = u1;
			u1 = u2;
			u2 = u3;
			u3 = u_new;
			cout << i << ":";
			for (int j = 0; j < 6; j++)
			{
				cout << u3[j] << " ";
			}
			cout << endl;
		}
		return u3;
	}

}
