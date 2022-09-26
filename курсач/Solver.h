#pragma once
#include "MyMatrix.h"
class Solver
{
private:
    MyMatrix A;
    MyVector p,
        z,
        r,
        x0,
        Ar,
        y;
    std::vector<double> L, D, U;
    int N;
    int maxIter;
    double eps;
    int iter;
    double normR;
    double normB;
public:

    Solver(int size);
    Solver(std::string filename);
    Solver(MyMatrix _A);

    void CGM_LU(); // conjugate gradient method whith LU factorization
    void LOS_LU(); // locally optimal scheme whith LU factorization
    void FactLU(std::vector<double>& L, std::vector<double>& U, std::vector<double>& D);
    void Direct(std::vector<double>& L, std::vector<double>& D, MyVector& y, MyVector& b);
    void Direct(std::vector<double>& L, MyVector& y, MyVector& b);
    void Reverse(std::vector<double>& U, MyVector& x, MyVector& y);
    void Reverse(std::vector<double>& U, std::vector<double>& D, MyVector& x, MyVector& y);
    void output(std::string filename);
    void getx0(std::vector<double>& x);

};
