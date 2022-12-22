#pragma once
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include <unordered_set>
#include <sstream>
#include<functional>

using namespace std;

//typedef // создаем новый прототип (в данном случае указатель на функцию)
//double // возвращаемое значение (такое же как в функциях)
//(*func_2coord_to_1) // имя прототипа (в коде употребляется без звездочки)
//(double, double); // список параметров (такое же как в функциях)

typedef std::function<double(double ksi, double eta)> func_2coord_to_1;

struct Node2D
{
	double x;
	double y;

	Node2D() { x = 0; y = 0; }
	Node2D(double _x, double _y) { x = _x; y = _y; }
};

//#include "includes.h"
//#include "FEM.h"

/*struct sum_elements3D  // последовательность заранее посчитанных psi, таких что x = SUM (x_i * psi_i)
{
	vector<double> psi;
	sum_elements3D() { psi.resize(8); }
};*/

struct Gauss_point3D
{
	double coord[3];
	double koef;
	double derivative[8][3]; // d_ksi, d_eta, d_teta  для этого coord[3]
	    // derivative_basis[mu(node_id)](ksi_x)* basis[nu(node_id)](ksi_y)* basis[teta(node_id)](ksi_z)
		// basis[mu(node_id)](ksi_x)* derivative_basis[nu(node_id)](ksi_y)* basis[teta(node_id)](ksi_z)
		// basis[mu(node_id)](ksi_x)* basis[nu(node_id)](ksi_y)* derivative_basis[teta(node_id)](ksi_z)
	double basis_func[8]; // basis[mu(node_id)](ksi_x)* basis[nu(node_id)](ksi_y)* basis[teta(node_id)](ksi_z)
};

struct Gauss_point2D
{
	double coord[2];
	double koef;
	double derivative[4][2]; // d_ksi, d_eta, d_teta  для этого coord[3]
		// derivative_basis[mu(node_id)](ksi_x)* basis[nu(node_id)](ksi_y)* basis[teta(node_id)](ksi_z)
		// basis[mu(node_id)](ksi_x)* derivative_basis[nu(node_id)](ksi_y)* basis[teta(node_id)](ksi_z)
		// basis[mu(node_id)](ksi_x)* basis[nu(node_id)](ksi_y)* derivative_basis[teta(node_id)](ksi_z)
	double basis_func[4]; // basis[mu(node_id)](ksi_x)* basis[nu(node_id)](ksi_y)* basis[teta(node_id)](ksi_z)
};

//class Integrate_Gauss3Method3D
//{
//	double basic_coef[3] = { 5. / 9., 8. / 9., 5. / 9. };
//	double basic_point[3] = { 0.77459666924148337703585307995648, 0, -0.77459666924148337703585307995648 };
//	double CalcdetJ(int p_id, vector<Node3D> &element);
//	double CalcdetJ_face(int p_id, vector<Node3D>& face, int face_id);
//	double CalcJ_1(vector<vector<double>>& J_1, int p_id, vector<Node3D> &element);
//public:
//	Gauss_point3D intagrate_points_3D[3*3*3];
//	Gauss_point3D intagrate_points_face[6][3*3]; // intagrate_points_face[0] - для интегрирования по 0ой грани
//
//	int GetLocG(vector<Node3D> &element, vector<vector<double>>& locG); // без домножения на лямбда
//	int GetLocM(vector<Node3D>& element, vector<vector<double>>& locM); // без домножения на гамма
//	int GetLocM_face(vector<Node3D>& face, vector<vector<double>>& locM, int fc_id); // без домножения на гамма
//
//	void Init(vector<basis_func>& basis, vector<basis_func>& derivative_basis);
//
//	double IntegrateFunc(Node3D& from, Node3D& to, func_3coord_to_1& func);
//	//Integrate_Gauss3Method3D(vector<basis_func> &basis, vector<basis_func> &derivative_basis);
//};

class Integrate_Gauss3Method2D
{
	//vector<vector<sum_elements3D>> J_base;
	vector<double> phi_in_points, derivative_phi_in_points; // значения функций и производных в точках
	double basic_coef[3] = { 5. / 9., 8. / 9., 5. / 9. };
	double basic_point[3] = { 0.77459666924148337703585307995648, 0, -0.77459666924148337703585307995648 };
	//double CalcdetJ(int p_id, vector<Node2D>& element);
	//double CalcJ_1(vector<vector<double>>& J_1, int p_id, vector<Node2D>& element);
public:
	vector<Gauss_point2D> intagrate_points_2D;

	//int GetLocG(vector<Node2D>& element, vector<vector<double>>& locG); // без домножения на лямбда
	//int GetLocM(vector<Node2D>& element, vector<vector<double>>& locM); // без домножения на гамма
	//void Init(vector<basis_func>& basis, vector<basis_func>& derivative_basis);

	void Init();
	double IntegrateFunc(Node2D& from, Node2D& to, func_2coord_to_1& func);

	//Integrate_Gauss3Method3D(vector<basis_func> &basis, vector<basis_func> &derivative_basis);
};