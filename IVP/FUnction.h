#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

template <typename T>
vector<T> operator*(const T& c, const vector<T>& u)
{
	vector<T> uu;
	for (int i = 0; i < u.size(); i++)
	{
		uu.push_back(c * u[i]);
	}
	return uu;
}

template <typename T>
vector<T> operator+(const vector<T>& u1, const vector<T>& u2)
{
	vector<T> uu;
	for (int i = 0; i < u1.size(); i++)
	{
		uu.push_back(u1[i] + u2[i]);
	}
	return uu;
}

template <typename T>
vector <T> operator-(const vector<T>& u1, const vector<T>& u2)
{
	vector<T> uu;
	for (int i = 0; i < u1.size(); i++)
	{
		uu.push_back(u1[i] - u2[i]);
	}
	return uu;
}

vector<double> f(const vector<double>& u, double t)
{
	vector<double> v;
	double mu = 0.01227747;
	v.push_back(u[3]);
	v.push_back(u[4]);
	v.push_back(u[5]);
	v.push_back(2 * u[4] + u[0] - (mu * (u[0] + mu - 1)) / pow((pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2)), 1.5) - ((1 - mu) * (u[0] + mu)) / pow((pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2)), 1.5));
	v.push_back(-2 * u[3] + u[1] - (mu * u[1]) / pow((pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2)), 1.5) - ((1 - mu) * u[1]) / pow((pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2)), 1.5));
	v.push_back(-(mu * u[2]) / pow((pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu - 1, 2)), 1.5) - ((1 - mu) * u[2]) / pow((pow(u[1], 2) + pow(u[2], 2) + pow(u[0] + mu, 2)), 1.5));
	/*for (int i = 0; i < v.size(); i++)
	{
		cout << v[i] << endl;
	}*/
	return v;
}
