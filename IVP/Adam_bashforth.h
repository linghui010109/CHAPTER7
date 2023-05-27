#include <iostream>
#include <cmath>
#include <vector>
#include "Classical_RK.h"
//#include "Function.h"

using namespace std;

 class Adam_Bashforth
 {
 private:
	 int n;
	 double T, mu;
 public:
	 vector<double> Adam_Bashforth_solver(int n,vector<double> &u,double T,int p);
 };

 vector<double> Adam_Bashforth::Adam_Bashforth_solver(int n,vector<double> &u, double T,int p)
 {
	 classical_RK r1;
	 double k = 1.0 * T / n;
	 
	 if (p == 1)
	 {
		 vector<double> u_new;
		 for (int i = 1; i < n; i++)
		 {
			 u_new = u + k * f(u, k * i);
			 u = u_new;
			 cout << i << ":";
			 for (int j = 0; j < 6; j++)
			 {
				 cout << u[j] << " ";
			 }
			 cout << endl;
		 }
		 return u;
	 }

	 if (p == 2)
	 {
		 vector<double> uu=r1.classical_RK_solver(n, u, T,1);
		 vector<double> u_new;
		 for (int i = 2; i < n; i++)
		 {
			 u_new = uu + k * (-0.5 * f(u, 0) + 1.5 * f(uu, 0));
			 u = uu;
			 uu = u_new;
			 cout << i << ":";
			 for (int j = 0; j < 6; j++)
			 {
				 cout << uu[j] << " ";
			 }
			 cout << endl;
		 }
		 return uu;
	 }
	 if (p == 3)
	 {
		 vector<double> u1 = r1.classical_RK_solver(n, u, T,1);
		 vector<double> u2 = r1.classical_RK_solver(n, u, T,2);
		 vector<double> u_new;
		 for (int i = 3; i < n; i++)
		 {
			 u_new = u2 + k * (5.0 / 12.0 * f(u, 0) - 4.0 / 3.0 * f(u1, 0) + 23.0 / 12.0 * f(u2, 0));
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

	 if (p == 4)
	 {
		 vector<double> u1 = r1.classical_RK_solver(n, u, T,1);
		 vector<double> u2 = r1.classical_RK_solver(n, u, T,2);
		 vector<double> u3 = r1.classical_RK_solver(n, u, T,3);
		 vector<double> u_new;
		 for (int i = 4; i < n; i++)
		 {
			 u_new = u3 + k * (-9.0 / 24.0 * f(u, 0) + 37.0 / 24.0 * f(u1, 0) - 59.0 / 24.0 * f(u2, 0) + 55.0 / 24.0 * f(u3, 0));
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
