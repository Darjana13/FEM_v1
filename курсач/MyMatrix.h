#pragma once
#include "MyVector.h"

class MyMatrix
{
public:
    std::vector<double> di;
    std::vector<double> al;
    std::vector<double> au;
    MyVector b;
    std::vector<int> ja;
    std::vector<int> ia;
    int N;

    MyMatrix(void);
    void ReadMatrix(int size);
    void Ax(std::vector<double>& X, std::vector<double>& Y); // ��������� ������� �� ������ y = Ax
    void Ax(MyVector& X, MyVector& Y); // ��������� ������� �� ������ y = Ax
    void ATx(MyVector& X, MyVector& Y); // ��������� ����������������� ������� �� ������ y = A^T*x
    MyMatrix& operator + (MyMatrix B);
    MyMatrix operator * (const double a);
    MyMatrix& operator=(const MyMatrix& B);

    void Show();

};