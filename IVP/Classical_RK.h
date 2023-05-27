#include <iostream>
#include <cmath>
#include <vector>
#include "Function.h"

using namespace std;

class classical_RK
{
private:
	//int n;
	//double T;
public:
	vector<double> classical_RK_solver(int n, vector<double> u, double T, int p);
};

vector<double> classical_RK::classical_RK_solver(int n, vector<double> u, double T, int p)
{
	double k = 1.0 * T / n;
	//vector<double> u = v;
	vector<double> u_new;
	for (int i = 0; i < p; i++)
	{
		vector<double> u_1, u_2, u_3, u_4;
		u_1 = f(u, 0);
		u_2 = f(u + 0.5 * k * u_1, 0.5 * k);
		u_3 = f(u + 0.5 * k * u_2, 0.5 * k);
		u_4 = f(u + k * u_3, k);
		u_new = u + ((k / 6) * (u_1 + 2.0 * u_2 + 2.0 * u_3 + u_4));
		u = u_new;
	}
	/*for (int i = 0; i < u.size(); i++)
	{
		cout << u[i] << endl;
	}*/
	return u;
}

/*vector<double> classical_RK::classical_RK_solver(int n, vector<double> u, double T)
{
	double k = 1.0 * T / n;
	//vector<double> u = v;
	vector<double> u_new;
	for (int i = 0; i < n; i++)
	{
		vector<double> u_1, u_2, u_3, u_4;
		u_1 = f(u, 0);
		u_2 = f(u + 0.5 * k * u_1, 0.5 * k);
		u_3 = f(u + 0.5 * k * u_2, 0.5 * k);
		u_4 = f(u + k * u_3, k);
		u_new = u + ((k / 6) * (u_1 + 2.0 * u_2 + 2.0 * u_3 + u_4));
		u = u_new;
	}
	/*for (int i = 0; i < u.size(); i++)
	{
		cout << u[i] << endl;
	}
	return u;
}*/
