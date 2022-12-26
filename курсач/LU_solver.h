#pragma once
#include <vector>
#include <string>
#include <fstream>

using namespace std;

class LU_solver
{
	const double eps = 10E-20;

public:
	int N = 0; // ������ �������
	vector<double> di; // ���������
	vector<int> ia; // ������� �������
	vector<double> au; // ������� �����������
	vector<double> al; // ������ �����������
	vector<double> b; // ������

	vector<double> di_LU; // ���������
	vector<double> au_LU; // ������� �����������
	vector<double> al_LU; // ������ �����������
	vector<double> y; // ������


	void Init(vector<int>& _ia, vector<int>& _ja, vector<double>& _au, vector<double>& _al, vector<double>& _di, vector<double>& _b);
	void Read(string path);
	void revers();
	void forward_stroke();
	bool LU();
	void Recalc(vector<double>& _b, vector<double>& _x);
	void Calc(vector<double>& _x);

};
