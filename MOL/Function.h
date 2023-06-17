#include<iostream>
#include<cmath>
#include<vector>
#include<cstring>
#include<fstream>
#include"Eigen/Dense"

using namespace std;

const double _epsL = 10 * std::numeric_limits<double>::epsilon();

class Function {
public:
    virtual double operator () (const double& x) const = 0;
};

class Func1 : public Function {
public:
    double operator () (const double& x) const {
        return exp(-20.0 * pow(x - 2.0, 2)) + exp(-pow(x - 5.0, 2));
    }
} func1;

#pragma once
