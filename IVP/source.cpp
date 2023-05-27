#include <iostream>
#include <cmath>
#include <vector>
//#include "Adam_bashforth.h"
//#include "Newton.h"
//#include "Function.h"
#include "BDFs.h"


using namespace std;

int main()
{
	
	vector<double> u;
    vector<double> initialvalue;
	double T = 17.06521656015796;
	double uu[6] = { 0.994,0,0,0,-2.0015851063790825224,0 };
	for (int i = 0; i < 6; i++)
	{
        initialvalue.push_back(0);
		u.push_back(uu[i]);
        //cout << u[i] << endl;
	}
    double k = T/2;

    //classical_RK r1;
	//r1.classical_RK_solver(10, u, T,10);

    //Adam_Bashforth a1;//p=1,2,3,4
    //a1.Adam_Bashforth_solver(10, u, T, 4);
    
   
    //Adam_Moulton m1;//p=2,3,4,5
    //m1.Adam_Moulton_solver(10, u, T, 5);

    BDFs B1;//p=1,2,3,4
    B1.BDFs_solver(10, u, T, 4);

	return 0;
}
/*
int main()
{
    vector<vector<int>> nums1 = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9} };
    vector<vector<int>> nums2 = {
        {1, 1, 1},
        {2, 2, 2},
        {3, 3, 3} };
    vector<vector<int>> nums3 = nums1 + nums2;

    for (int i = 0; i < nums3.size(); i++)
    {
        for (int j = 0; j < nums3[0].size(); j++)
        {
            cout << nums3[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}*/
