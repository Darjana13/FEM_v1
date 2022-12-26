#pragma once
#include <vector>
#include <string>
#include <fstream>

using namespace std;

class LU_solver
{
	const double eps = 10E-20;

public:
	int N = 0; // Размер матрицы
	vector<double> di; // Диагональ
	vector<int> ia; // Портрет матрицы
	vector<double> au; // верхний треугольник
	vector<double> al; // нижний треугольник
	vector<double> b; // Вектор

	vector<double> di_LU; // Диагональ
	vector<double> au_LU; // верхний треугольник
	vector<double> al_LU; // нижний треугольник
	vector<double> y; // Вектор


	void Init(vector<int>& _ia, vector<int>& _ja, vector<double>& _au, vector<double>& _al, vector<double>& _di, vector<double>& _b);
	void Read(string path);
	void revers();
	void forward_stroke();
	bool LU();
	void Recalc(vector<double>& _b, vector<double>& _x);
	void Calc(vector<double>& _x);

};
