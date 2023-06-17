#include"Header.h"
#include"Function.h"

int main()
{
//	A a1;
	//a1.create_A(5, 1, 1, 1);
	double h = 0.05;
	int n = 20;
	int T = 17;
	double k = 0.8 * h;
	double a = 1.0;
	Lax_Friedrichs l1;
	Func1 f1;
	MatrixXd u_0(n, 1);
	for (int i = 0; i < n; i++)
	{
		u_0(i, 0) = f1((i + 1) * h*10.0);
		//cout << f1(2.0) << endl;
		cout << (i + 1) * h * 10.0 <<" " << u_0(i, 0) << endl;
	}
	
	//MatrixXd u=l1.solve_lax_Friedrichs(n, a, h, k, T, u_0);
	//cout << u << endl;

	//Lax_Wendroff l2;
	//MatrixXd u2 = l2.solve_lax_Wendroff(n, a, h, k, T, u_0);
	//cout << u2 << endl;

	upwind l3;
	MatrixXd u3 = l3.solve_upwind(n, a, h, k, T, u_0);
	cout << u3 << endl;
	return 0;

}
