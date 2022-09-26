#pragma once
#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
class MyVector
{
public:
    std::vector<double> vect;
    MyVector();
    void Size(int N);
    void ReadVector(std::string filename);
    MyVector& operator=(const MyVector& a);
    double operator* (const MyVector& a);
    MyVector& operator* (const double a);
    MyVector& operator-(MyVector& a);
    MyVector& operator+(MyVector& a);
    double Norm();
};


