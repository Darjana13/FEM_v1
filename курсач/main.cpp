#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include "Solver.h"
#include "Generate.h"
#include "Gauss.h"
#include <Windows.h>
#include "LU_solver.h"

void ClearFolder(string path)
{
	string cmd = "del /s /q " + path + "*.*";
	system(cmd.c_str());
}

template <typename T>
string NumberToString(T Number)
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}

typedef // создаем новый прототип (в данном случае указатель на функцию)
double // возвращаемое значение (такое же как в функциях)
(*basis_func) // имя прототипа (в коде употребляется без звездочки)
(double, double, double); // список параметров (такое же как в функциях)

double X1(double x0, double x1, double x)
{
	return (x1 - x) / (x1 - x0);
}
double X2(double x0, double x1, double x)
{
	return (x - x0) / (x1 - x0);
}

double dX1(double x0, double x1, double x)
{
	return (-1.0) / (x1 - x0);
}
double dX2(double x0, double x1, double x)
{
	return (1.0) / (x1 - x0);
}

vector<basis_func> basis_1D = { X1, X2 }, dbasis_1D = { dX1, dX2 };

MyVector q1, q2, q3, q4;
Integrate_Gauss3Method2D gauss;

struct node
{
	double r;
	double z;
};

struct material
{
	double lambda;
	int gamma_id; // number of gamma function
};

struct element
{
	std::vector<int> node_loc;
	int mater;
	int f_id;
};

std::vector<node> all_nodes;         // все узлы в порядке глобальной нумерации
std::vector<element> all_elems;      // все элементы в прядке глобальной нумерации
std::vector<material> all_materials; // все материалы по индексам
std::vector<std::pair<int, std::vector<int>>> S1;      // S1[i][j] на j-ом узле заданы краевые 1 рода
std::vector<std::pair<int, std::vector<int>>> S2_r;    // граница параллельна r на j-ом узле заданы краевые 2 рода
std::vector<std::pair<int, std::vector<int>>> S2_z;    // граница параллельна z на j-ом узле заданы краевые 2 рода
std::vector<std::pair<int, std::vector<int>>> S3_r;    // граница параллельна r на j-ом узле заданы краевые 3 рода
std::vector<std::pair<int, std::vector<int>>> S3_z;    // граница параллельна z на j-ом узле заданы краевые 3 рода

std::vector<double> time_grid; // сетка времени
int i_t = 0; // текущий временной слой

double gamma(double r, double z, int gam_id) // значение гамма по индексу gam_id 
{
	switch (gam_id)
	{
	case 0:
		return 
			4200000
			;
	case 1:
		return 
			2430000
			;
	default:
		std::cout << "can't find gamma № " << gam_id << "\n";
		break;
	}
}
double z_max_abs = 0.15, r_max_abs = 0.1; // !!!

//int GetVectorV(std::vector<double>& v, double r, double z, double z_max, double r_max)
//{
//	double R = r_max, H = z_max, v_max = 0.02;
//	//vector <double> v;
//
//	// 1___________________________
//	if (z_max - r >= z && z >= r && r <= 0.5 * r_max)
//	{
//		//cout « "1" « "\n";
//		if (r < R / 4)
//		{
//			v[0] = 0;
//			v[1] = -v_max * r / (R / 4);
//		}
//		else
//		{
//			v[0] = 0;
//			v[1] = -v_max * (R / 2 - r) / (R / 4);
//		}
//		//cout « v[0] « " " « v[1];
//	}
//	// 2___________________________
//	if (z_max - r < z && z > r && r >= 0.5 * z_max) //(z_max - r < z && z > r && z <= 0.5 * z_max)
//	{
//		//cout « "2" « "\n";
//		if (r > 3 * H / 4)
//		{
//			v[1] = 0;
//			v[0] = -v_max * (H - z) / (H / 4);
//		}
//		else
//		{
//			v[1] = 0;
//			v[0] = -v_max * (z - H / 2) / (H / 4);
//		}
//		//cout « v[0] « " " « v[1];
//	}
//	//3 _______________
//	if (z_max - r <= z && z <= r && r >= 0.5 * r_max)
//
//	{
//		//cout « "3" « "\n";
//		if (r > 3 * R / 4)
//		{
//			v[0] = 0;
//			v[1] = v_max * (R - r) / (R / 4);
//		}
//		else
//		{
//			v[0] = 0;
//			v[1] = v_max * (r - R / 2) / (R / 4);
//		}
//		//cout « v[0] « " " « v[1];
//	}
//	//4 ________________
//	if (z_max - r <= r && z < r && z <= 0.5 * z_max)
//	{
//		//cout « "4" « "\n";
//		if (z < H / 4)
//		{
//			v[1] = 0;
//			v[0] = v_max * r / (H / 4);
//		}
//		else
//		{
//			v[1] = 0;
//			v[0] = v_max * (H / 2 - z) / (H / 4);
//		}
//		//cout « v[0]« " " « v[1];
//	}
//
//	return 0;
//}

int GetVectorV(std::vector<double>& v, double r, double z, double z_max, double r_max)
{
	//v[0] = 0;
	//v[1] = 0;
	//return -1;

	double R = r_max, H = z_max, v_max = 0.001;
	z -= 0.01;
	if (r > r_max || z < 0)
	{
		v[0] = 0;
		v[1] = 0;

		return -1;
	}
	//vector <double> v;
	int num = 0;
	// 1___________________________
	if (z - 1.5*r >= 0 && z + 1.5*r - z_max < 0)
	{
		num = 1;
		//cout << "1";
		if (r < R / 4)
		{
			v[0] = 0;
			v[1] = -v_max * r / (R / 4);
		}
		else
		{
			v[0] = 0;
			v[1] = -v_max * (R / 2 - r) / (R / 4);
		}
		//cout << v[0] << " " << v[1];
	}
	// 2___________________________
	if (z - 1.5 * r >= 0 && z + 1.5 * r - z_max >= 0) //(z_max - r < z && z > r && z <= 0.5 * z_max)
	{
		num = 2;
		//cout << "2";
		if (z > 3 * H / 4)
		{
			v[1] = 0;
			v[0] = -v_max * (H - z) / (H / 4);
		}
		else
		{
			v[1] = 0;
			v[0] = -v_max * (z - H / 2) / (H / 4);
		}
		//cout << v[0] << " " << v[1] << '\n';
	}
	//3 _______________
	if (z - 1.5 * r < 0 && z + 1.5 * r - z_max >= 0)

	{
		num = 3;

		//cout << "3";
		if (r > 3 * R / 4)
		{
			v[0] = 0;
			v[1] = v_max * (R - r) / (R / 4);
		}
		else
		{
			v[0] = 0;
			v[1] = v_max * (r - R / 2) / (R / 4);
		}
		if (v[1] < 0)
			cout << v[0] << " " << v[1] << " r " << r << " z " << z << '\n';
		//cout << v[0] << " " << v[1] << '\n';
	}
	//4 ________________
	if (z - 1.5 * r < 0 && z + 1.5 * r - z_max < 0)
	{
		num = 4;

		//cout << "4" << "\n";
		if (z < H / 4)
		{
			v[1] = 0;
			v[0] = v_max * z / (H / 4);
		}
		else
		{
			v[1] = 0;
			v[0] = v_max * (H / 2 - z) / (H / 4);
		}
		//cout << v[0] << " " << v[1] << '\n';
	}

	return num;
}
double beta(double r, double z, int beta_id) // считаем, что бета везде одинаковая
{
	switch (beta_id)
	{
	case 0:
		return 1;
	default:
		std::cout << "can't find gamma № " << beta_id << "\n";
		break;
	}
}



double func_f(double r, double z, int f_id) // значение f по индексу f_id 
{
	double t = time_grid[i_t];
	vector<double> v(2);
	GetVectorV(v, r, z, z_max_abs, r_max_abs);
	switch (f_id)
	{
	case 0:
		return 0
			//-0.56 / r + 4200000 * (1 + v[0] * r - v[1]);
			;
	default:
		std::cout << "can't find f № " << f_id << "\n";
		break;
	}
}

double func_S(double r, double z, int s_id) // значение краевого S по индексу f_id
{
	double t = time_grid[i_t];
	switch (s_id)
	{
	case 0:
		return //r * r + 3 * z - 5 * t
			52613.0 / 209.0 // 2000 / (pi * 0.11 * 0.11) / 209.0
			;
	case 1:// 2_z
		return 
			//2 * r * r
			r + z + t
			;
	case 2: // 2_r
		return 3
			;
	case 3: // 3_z
		return 1 / beta(r, z, 0) * 2 * r * r + r * r + 3 * z - 5 * t
			;
	case 4: // 3_r
		return 1 / beta(r, z, 0) * 3 + r * r + 3 * z - 5 * t
			;
	default:
		std::cout << "can't find S № " << s_id << "\n";
		break;
	}
}

int Input() // чтение данных
{
	int N, Nmat, Kel, NS1, Ntime, NS;
	std::ifstream in;

	in.open("info.txt");
	//in >> N >> Nmat >> Kel /* >> NS1*/;
	in.close();

	in.open("rz.txt");
	in >> N;
	all_nodes.resize(N);
	for (int i = 0; i < N; i++)
	{
		in >> all_nodes[i].r >> all_nodes[i].z;
	}
	in.close();

	in.open("time.txt");
	in >> Ntime;
	time_grid.resize(Ntime);
	for (int i = 0; i < Ntime; i++)
	{
		in >> time_grid[i];
	}
	in.close();

	// для тестирования
	/*std::ofstream temp("q0.txt");
	temp.precision(15);
	for (int k = 0; k < 1; k++)
	{
		i_t = k;
		for (int i = 0; i < N; i++)
		{
			//temp << func_S(all_nodes[i].r, all_nodes[i].z, 0) << " ";
			temp << 20.0 << " ";
		}
		temp << "\n";
	}
	temp.close();
	// конец для тестирования
	in.open("q0.txt");
	q1.Size(N);
	for (int i = 0; i < N; i++)
	{
		in >> q1.vect[i];
	}
	in.close();*/


	// для тестирования
	std::ofstream temp("q0 q1 q2.txt");
	temp.precision(15);
	for (int k = 0; k < 3; k++)
	{
	    i_t = k;
	    for (int i = 0; i < N; i++)
	    {
	        //temp << func_S(all_nodes[i].r, all_nodes[i].z, 0) << " ";
			temp << 20.0 << " ";

	    }
	    temp << "\n";
	}
	temp.close();
	// конец для тестирования

	in.open("q0 q1 q2.txt");
	q1.Size(N);
	for (int i = 0; i < N; i++)
	{
	    in >> q1.vect[i];
	}
	q2.Size(N);
	for (int i = 0; i < N; i++)
	{
	    in >> q2.vect[i];
	}
	q3.Size(N);
	for (int i = 0; i < N; i++)
	{
	    in >> q3.vect[i];
	}
	q4.Size(N);
	in.close();

	in.open("S1.txt");
	in >> NS1;
	S1.resize(NS1);
	for (int i = 0; i < NS1; i++)
	{
		int size;
		in >> size >> S1[i].first;
		S1[i].second.resize(size);
		for (int j = 0; j < size; j++)
		{
			in >> S1[i].second[j];
		}
	}
	in.close();

	in.open("S2_r.txt");
	in >> NS;
	S2_r.resize(NS);
	for (int i = 0; i < NS; i++)
	{
		int size;
		in >> size >> S2_r[i].first;
		S2_r[i].second.resize(size);
		for (int j = 0; j < size; j++)
		{
			in >> S2_r[i].second[j];
		}
	}
	in.close();

	in.open("S2_z.txt");
	in >> NS;
	S2_z.resize(NS);
	for (int i = 0; i < NS; i++)
	{
		int size;
		in >> size >> S2_z[i].first;
		S2_z[i].second.resize(size);
		for (int j = 0; j < size; j++)
		{
			in >> S2_z[i].second[j];
		}
	}
	in.close();

	in.open("S3_r.txt");
	in >> NS;
	S3_r.resize(NS);
	for (int i = 0; i < NS; i++)
	{
		int size;
		in >> size >> S3_r[i].first;
		S3_r[i].second.resize(size);
		for (int j = 0; j < size; j++)
		{
			in >> S3_r[i].second[j];
		}
	}
	in.close();

	in.open("S3_z.txt");
	in >> NS;
	S3_z.resize(NS);
	for (int i = 0; i < NS; i++)
	{
		int size;
		in >> size >> S3_z[i].first;
		S3_z[i].second.resize(size);
		for (int j = 0; j < size; j++)
		{
			in >> S3_z[i].second[j];
		}
	}
	in.close();

	in.open("material.txt");
	in >> Nmat;
	all_materials.resize(Nmat);
	for (int i = 0; i < Nmat; i++)
	{
		in >> all_materials[i].lambda >> all_materials[i].gamma_id;
	}
	in.close();

	in.open("elem.txt");
	in >> Kel;
	all_elems.resize(Kel);
	for (int i = 0; i < Kel; i++)
	{
		all_elems[i].node_loc.resize(4);
		in >> all_elems[i].node_loc[0] >> all_elems[i].node_loc[1]
			>> all_elems[i].node_loc[2] >> all_elems[i].node_loc[3]
			>> all_elems[i].mater >> all_elems[i].f_id;
	}
	in.close();

	return 0;
}
double GetG_Loc(double rp, double lambda, double hr, double hz,
	std::vector<std::vector<double>>& G_loc) // получение локальной G
{
	double a1 = (lambda * hz * rp) / (6 * hr),
		a2 = (lambda * hz) / (12),
		a3 = (lambda * hr * rp) / (6 * hz),
		a4 = (lambda * hr * hr) / (12 * hz);
	G_loc[0][0] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
	G_loc[0][1] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
	G_loc[0][2] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
	G_loc[0][3] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;

	G_loc[1][0] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
	G_loc[1][1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
	G_loc[1][2] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
	G_loc[1][3] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;

	G_loc[2][0] = 1 * a1 + 1 * a2 - 2 * a3 - 1 * a4;
	G_loc[2][1] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
	G_loc[2][2] = 2 * a1 + 2 * a2 + 2 * a3 + 1 * a4;
	G_loc[2][3] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;

	G_loc[3][0] = -1 * a1 - 1 * a2 - 1 * a3 - 1 * a4;
	G_loc[3][1] = 1 * a1 + 1 * a2 - 2 * a3 - 3 * a4;
	G_loc[3][2] = -2 * a1 - 2 * a2 + 1 * a3 + 1 * a4;
	G_loc[3][3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
	return 0;
}
double GetM_Loc(double rp, double zs, int gam, double hr, double hz,
	std::vector<std::vector<double>>& M_loc) // прибавление локальной М
{
	double g1 = gamma(rp, zs, gam),
		g2 = gamma(rp + hr, zs, gam),
		g3 = gamma(rp, zs + hz, gam),
		g4 = gamma(rp + hr, zs + hz, gam);
	M_loc[0][0] = hr * (
		g1 * (rp / 4 + hr / 20) * hz / 4 +
		g2 * (rp / 12 + hr / 30) * hz / 4 +
		g3 * (rp / 4 + hr / 20) * hz / 12 +
		g4 * (rp / 12 + hr / 30) * hz / 12);
	M_loc[0][1] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 4 +
		g2 * (rp / 12 + hr / 20) * hz / 4 +
		g3 * (rp / 12 + hr / 30) * hz / 12 +
		g4 * (rp / 12 + hr / 20) * hz / 12);
	M_loc[0][2] = hr * (
		g1 * (rp / 4 + hr / 20) * hz / 12 +
		g2 * (rp / 12 + hr / 30) * hz / 12 +
		g3 * (rp / 4 + hr / 20) * hz / 12 +
		g4 * (rp / 12 + hr / 30) * hz / 12);
	M_loc[0][3] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 12 +
		g2 * (rp / 12 + hr / 20) * hz / 12 +
		g3 * (rp / 12 + hr / 30) * hz / 12 +
		g4 * (rp / 12 + hr / 20) * hz / 12);
	M_loc[1][0] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 4 +
		g2 * (rp / 12 + hr / 20) * hz / 4 +
		g3 * (rp / 12 + hr / 30) * hz / 12 +
		g4 * (rp / 12 + hr / 20) * hz / 12);
	M_loc[1][1] = hr * (
		g1 * (rp / 12 + hr / 20) * hz / 4 +
		g2 * (rp / 4 + hr / 5) * hz / 4 +
		g3 * (rp / 12 + hr / 20) * hz / 12 +
		g4 * (rp / 4 + hr / 5) * hz / 12);
	M_loc[1][2] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 12 +
		g2 * (rp / 12 + hr / 20) * hz / 12 +
		g3 * (rp / 12 + hr / 30) * hz / 12 +
		g4 * (rp / 12 + hr / 20) * hz / 12);
	M_loc[1][3] = hr * (
		g1 * (rp / 12 + hr / 20) * hz / 12 +
		g2 * (rp / 4 + hr / 5) * hz / 12 +
		g3 * (rp / 12 + hr / 20) * hz / 12 +
		g4 * (rp / 4 + hr / 5) * hz / 12);
	M_loc[2][0] = hr * (
		g1 * (rp / 4 + hr / 20) * hz / 12 +
		g2 * (rp / 12 + hr / 30) * hz / 12 +
		g3 * (rp / 4 + hr / 20) * hz / 12 +
		g4 * (rp / 12 + hr / 30) * hz / 12);
	M_loc[2][1] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 12 +
		g2 * (rp / 12 + hr / 20) * hz / 12 +
		g3 * (rp / 12 + hr / 30) * hz / 12 +
		g4 * (rp / 12 + hr / 20) * hz / 12);
	M_loc[2][2] = hr * (
		g1 * (rp / 4 + hr / 20) * hz / 12 +
		g2 * (rp / 12 + hr / 30) * hz / 12 +
		g3 * (rp / 4 + hr / 20) * hz / 4 +
		g4 * (rp / 12 + hr / 30) * hz / 4);
	M_loc[2][3] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 12 +
		g2 * (rp / 12 + hr / 20) * hz / 12 +
		g3 * (rp / 12 + hr / 30) * hz / 4 +
		g4 * (rp / 12 + hr / 20) * hz / 4);
	M_loc[3][0] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 12 +
		g2 * (rp / 12 + hr / 20) * hz / 12 +
		g3 * (rp / 12 + hr / 30) * hz / 12 +
		g4 * (rp / 12 + hr / 20) * hz / 12);
	M_loc[3][1] = hr * (
		g1 * (rp / 12 + hr / 20) * hz / 12 +
		g2 * (rp / 4 + hr / 5) * hz / 12 +
		g3 * (rp / 12 + hr / 20) * hz / 12 +
		g4 * (rp / 4 + hr / 5) * hz / 12);
	M_loc[3][2] = hr * (
		g1 * (rp / 12 + hr / 30) * hz / 12 +
		g2 * (rp / 12 + hr / 20) * hz / 12 +
		g3 * (rp / 12 + hr / 30) * hz / 4 +
		g4 * (rp / 12 + hr / 20) * hz / 4);
	M_loc[3][3] = hr * (
		g1 * (rp / 12 + hr / 20) * hz / 12 +
		g2 * (rp / 4 + hr / 5) * hz / 12 +
		g3 * (rp / 12 + hr / 20) * hz / 4 +
		g4 * (rp / 4 + hr / 5) * hz / 4);
	return 0;
}

double GetMG_Loc(double rp, double zs, int gam, double hr, double hz,
	std::vector<std::vector<double>>& MG_loc) // прибавление локальной М
{
	/*double g1 = gamma(rp, zs, gam),
		g2 = gamma(rp + hr, zs, gam),
		g3 = gamma(rp, zs + hz, gam),
		g4 = gamma(rp + hr, zs + hz, gam);*/

	Node2D from(rp, zs), to(rp + hr, zs + hz);
	func_2coord_to_1 phi;
	std::vector<double> v(2);
	GetVectorV(v, rp + hr / 2.0, zs + hz / 2.0, z_max_abs, r_max_abs);
	double gr = gamma(rp, zs, gam) * v[0];
	double gz = gamma(rp, zs, gam) * v[1];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			/*f = [i, x, y, z](double ksi, double eta, double zeta)
			{
				double J_1[3][3], det_J, xcrd(0.0), ycrd(0.0), zcrd(0.0);
				Jacobian3D(ksi, eta, zeta, x, y, z, J_1, det_J);
				for (int k = 0; k < 8; k++)
				{
					xcrd += x[k] * phi3d(k, ksi, eta, zeta);
					ycrd += y[k] * phi3d(k, ksi, eta, zeta);
					zcrd += z[k] * phi3d(k, ksi, eta, zeta);
				}
				return phi3d(i, ksi, eta, zeta) * func(xcrd, ycrd, zcrd) * abs(det_J);
			};*/
			phi = [i, j, from, to, gr, gz](double x, double y)
			{
				return (gr * dbasis_1D[j % 2](from.x, to.x, x) * basis_1D[j / 2](from.y, to.y, y) +
					gz * basis_1D[j % 2](from.x, to.x, x) * dbasis_1D[j / 2](from.y, to.y, y)) *
					basis_1D[i % 2](from.x, to.x, x) * basis_1D[i / 2](from.y, to.y, y) * x;
			};
			MG_loc[i][j] = gauss.IntegrateFunc(from, to, phi);
			//cout << MG_loc[i][j] << "\t";
		}
		//cout << "\n";
	}

	return 0;
}


int Getb_Loc(double rp, double zs, double hr, double hz,
	std::vector<double>& b_loc, int f_id) // получение локального b
{
	double f1 = func_f(rp, zs, f_id),
		f2 = func_f(rp + hr, zs, f_id),
		f3 = func_f(rp, zs + hz, f_id),
		f4 = func_f(rp + hr, zs + hz, f_id);
	b_loc[0] =
		f1 * (hr * hz / 3 * (rp / 3 + hr / 12)) +
		f2 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
		f3 * (hr * hz / 6 * (rp / 3 + hr / 12)) +
		f4 * (hr * hz / 6 * (rp / 6 + hr / 12));
	b_loc[1] =
		f1 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
		f2 * (hr * hz / 3 * (rp / 3 + hr / 4)) +
		f3 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
		f4 * (hr * hz / 6 * (rp / 3 + hr / 4));
	b_loc[2] =
		f1 * (hr * hz / 6 * (rp / 3 + hr / 12)) +
		f2 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
		f3 * (hr * hz / 3 * (rp / 3 + hr / 12)) +
		f4 * (hr * hz / 3 * (rp / 6 + hr / 12));
	b_loc[3] =
		f1 * (hr * hz / 6 * (rp / 6 + hr / 12)) +
		f2 * (hr * hz / 6 * (rp / 3 + hr / 4)) +
		f3 * (hr * hz / 3 * (rp / 6 + hr / 12)) +
		f4 * (hr * hz / 3 * (rp / 3 + hr / 4));
	return 0;
}

int Get_Loc(std::vector<std::vector<double>>& M_loc, std::vector<std::vector<double>>& G_loc, std::vector<std::vector<double>>& MG_loc,
	int el_id) // получение локальной матрицы А
{
	element el = all_elems[el_id];
	double hr = all_nodes[el.node_loc[1]].r - all_nodes[el.node_loc[0]].r,
		hz = all_nodes[el.node_loc[2]].z - all_nodes[el.node_loc[0]].z;

	// получить М
	GetM_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, all_materials[el.mater].gamma_id, hr, hz, M_loc);

	// получить G
	GetG_Loc(all_nodes[el.node_loc[0]].r, all_materials[el.mater].lambda, hr, hz, G_loc); // A_loc = G_loc 

	// получить МG, для конвекции
	GetMG_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, all_materials[el.mater].gamma_id, hr, hz, MG_loc);
	return 0;
}

int Get_Loc_b(std::vector<double>& b_loc,
	int el_id) // получить вектор правой части
{
	element el = all_elems[el_id];
	double hr = all_nodes[el.node_loc[1]].r - all_nodes[el.node_loc[0]].r,
		hz = all_nodes[el.node_loc[2]].z - all_nodes[el.node_loc[0]].z;
	// получить f
	Getb_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, hr, hz, b_loc, el.f_id);
	return 0;
}
int GeneratePortrait(MyMatrix& A,
	int N, int Kel) // генерация портрета
{
	std::vector<int>* ia = &A.ia,
		* ja = &A.ja;
	ia->resize(N + 1);
	ja->resize(16 * Kel);
	std::vector<int> temp_list1(16 * Kel),
		temp_list2(16 * Kel);
	std::vector<int> listbeg(N);
	int listsize = 0;
	for (int i = 0; i < N; i++)
	{
		listbeg[i] = 0;
	}
	for (int ielem = 0; ielem < Kel; ielem++)
	{
		for (int i = 0; i < 4; i++) // NumberOfUnknowns(ielem)?
		{
			int k = all_elems[ielem].node_loc[i];
			for (int j = i + 1; j < 4; j++)// NumberOfUnknowns(ielem)?
			{
				int ind1 = k;
				int ind2 = all_elems[ielem].node_loc[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				int iaddr = listbeg[ind2];
				if (iaddr == 0)
				{
					listsize++;
					listbeg[ind2] = listsize;
					temp_list1[listsize] = ind1;
					temp_list2[listsize] = 0;
				}
				else
				{
					while (temp_list1[iaddr] < ind1 && temp_list2[iaddr] > 0)
					{
						iaddr = temp_list2[iaddr];
					}
					if (temp_list1[iaddr] > ind1)
					{
						listsize++;
						temp_list1[listsize] = temp_list1[iaddr];
						temp_list2[listsize] = temp_list2[iaddr];
						temp_list1[iaddr] = ind1;
						temp_list2[iaddr] = listsize;
					}
					else if (temp_list1[iaddr] < ind1)
					{
						listsize++;
						temp_list2[iaddr] = listsize;
						temp_list1[listsize] = ind1;
						temp_list2[listsize] = 0;
					}
				}
			}
		}
	}

	(*ia)[0] = 0;
	for (int i = 0; i < N; i++)
	{
		(*ia)[i + 1] = (*ia)[i];
		int iaddr = listbeg[i];
		while (iaddr != 0)
		{
			(*ja)[(*ia)[i + 1]] = temp_list1[iaddr];
			(*ia)[i + 1]++;
			iaddr = temp_list2[iaddr];
		}
	}

	ja->resize((*ia)[N]);
	return 0;
}

int AddLocal(std::vector<int>& iaM, std::vector<int>& jaM, std::vector<double>& diM,
	std::vector<double>& alM, std::vector<double>& auM,
	std::vector<std::vector<double>>& M_loc,
	int el_id)
	// внесение локальных A, b  в глобальную СЛАУ
{
	std::vector<int> L = all_elems[el_id].node_loc;
	int n_loc = all_elems[el_id].node_loc.size(); // размерность локальной матрицы
	for (int i = 0; i < n_loc; i++)
	{
		diM[L[i]] += M_loc[i][i];
		//diG[L[i]] += G_loc[i][i];
	}

	for (int i = 0; i < 4; i++)
	{
		int temp = iaM[L[i]];
		for (int j = 0; j < i; j++)
		{
			for (int k = temp; k < iaM[L[i] + 1]; k++)
			{
				if (jaM[k] == L[j])
				{
					alM[k] += M_loc[i][j];
					auM[k] += M_loc[j][i];
					//alG[k] += G_loc[i][j];
					//auG[k] += G_loc[j][i];
					k++;
					break;
				}
			}
		}
	}

	return 0;
}

int AddLocal_b(std::vector<double>& b, std::vector<double>& b_loc, int el_id)
// внесение локальных b  в глобальную СЛАУ
{
	std::vector<int> L = all_elems[el_id].node_loc;
	int k = all_elems[el_id].node_loc.size(); // размерность локальной матрицы
	for (int i = 0; i < k; i++)
	{
		b[L[i]] += b_loc[i];
	}
	return 0;
}

int SetS1(std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& di,
	std::vector<double>& al, std::vector<double>& au,
	std::vector<double>& b) // учет первых краевых
{
	int NS1 = S1.size();
	for (int i = 0; i < NS1; i++)
	{
		int s1_id = S1[i].first;
		for (int j = 0; j < S1[i].second.size(); j++)
		{
			int node_id = S1[i].second[j];
			di[node_id] = 1;
			b[node_id] = func_S(all_nodes[node_id].r, all_nodes[node_id].z, s1_id);
			for (int k = ia[node_id]; k < ia[node_id + 1]; k++)
			{
				al[k] = 0;
			}
			for (int k = 0; k < ja.size(); k++)
			{
				if (ja[k] == node_id)
				{
					au[k] = 0;
				}
			}
		}
	}
	return 0;
}



double GetM_Loc_dim2_r(double rp, double hr,
	std::vector<std::vector<double>>& M_loc) // локальная матрица для S3_r
{

	M_loc[0][0] = hr / 6 * (2 * rp + hr / 2);
	M_loc[0][1] = hr / 6 * (rp + hr / 2);
	M_loc[1][0] = hr / 6 * (rp + hr / 2);
	M_loc[1][1] = hr / 6 * (2 * rp + 3 * hr / 2);

	return 0;
}
int Getb_Loc_dim2_r(double rp, double zs, double hr,
	std::vector<double>& b_loc, int f_id) // получение локального b для S параллельных r
{
	double f1 = func_S(rp, zs, f_id),
		f2 = func_S(rp + hr, zs, f_id);
	b_loc[0] = f1 * (hr / 6 * (2 * rp + hr / 2)) + f2 * (hr / 6 * (rp + hr / 2));
	b_loc[1] = f1 * (hr / 6 * (rp + hr / 2)) + f2 * (hr / 6 * (2 * rp + 3 * hr / 2));
	return 0;
}
double GetM_Loc_dim2_z(double zp, double hz,
	std::vector<std::vector<double>>& M_loc) // локальная матрица для S3_z
{
	M_loc[0][0] = hz / 3;
	M_loc[0][1] = hz / 6;
	M_loc[1][0] = hz / 6;
	M_loc[1][1] = hz / 3;

	return 0;
}
int Getb_Loc_dim2_z(double rp, double zs, double hz,
	std::vector<double>& b_loc, int s_id) // получение локального b для S параллельных z
{
	double f1 = func_S(rp, zs, s_id);
	double f2 = func_S(rp, zs + hz, s_id);
	b_loc[0] = f1 * (hz / 3) + f2 * (hz / 6);
	b_loc[1] = f1 * (hz / 6) + f2 * (hz / 3);
	return 0;
}
int AddLocal_dim2(std::vector<int>& iaM, std::vector<int>& jaM, std::vector<double>& diM,
	std::vector<double>& alM, std::vector<double>& auM,
	std::vector<std::vector<double>>& M_loc,
	int node1, int node2) // внесение локальной матрицы в глобальную для одномерного случая
{
	std::vector<int> L(2);
	L[0] = node1;
	L[1] = node2;
	int n_loc = 2; // размерность локальной матрицы
	for (int i = 0; i < n_loc; i++)
	{
		diM[L[i]] += M_loc[i][i];
	}

	for (int i = 0; i < 2; i++)
	{
		int temp = iaM[L[i]];
		for (int j = 0; j < i; j++)
		{
			for (int k = temp; k < iaM[L[i] + 1]; k++)
			{
				if (jaM[k] == L[j])
				{
					alM[k] += M_loc[i][j];
					auM[k] += M_loc[j][i];
					k++;
					break;
				}
			}
		}
	}

	return 0;
}
int Set_S2(MyMatrix& MS) // учет вторых краевых MS.b - вектор вклада от краевых
{
	std::vector<double> b_loc(2);
	int NS2 = S2_r.size();
	for (int i = 0; i < NS2; i++)
	{
		int s2_id = S2_r[i].first;
		for (int j = 0; j < S2_r[i].second.size() - 1; j++)
		{
			int node_id1 = S2_r[i].second[j],
				node_id2 = S2_r[i].second[j + 1];
			double hr = all_nodes[node_id2].r - all_nodes[node_id1].r;
			Getb_Loc_dim2_r(all_nodes[node_id1].r, all_nodes[node_id1].z, hr, b_loc, s2_id);
			MS.b.vect[node_id1] += b_loc[0];
			MS.b.vect[node_id2] += b_loc[1];
		}
	}
	NS2 = S2_z.size();
	for (int i = 0; i < NS2; i++)
	{
		int s2_id = S2_z[i].first;
		for (int j = 0; j < S2_z[i].second.size() - 1; j++)
		{
			int node_id1 = S2_z[i].second[j],
				node_id2 = S2_z[i].second[j + 1];
			double hz = all_nodes[node_id2].z - all_nodes[node_id1].z;
			Getb_Loc_dim2_z(all_nodes[node_id1].r, all_nodes[node_id1].z, hz, b_loc, s2_id);
			MS.b.vect[node_id1] += b_loc[0];
			MS.b.vect[node_id2] += b_loc[1];
		}
	}
	return 0;
}
int Set_S3(MyMatrix& MS, bool flag) // MS.b - вектор вклада от краевых
{
	std::vector<double> b_loc(2);
	std::vector<std::vector<double>> M_loc(2);
	M_loc[0].resize(2);
	M_loc[1].resize(2);

	int NS2 = S3_r.size();
	for (int i = 0; i < NS2; i++)
	{
		int s3_id = S3_r[i].first;
		for (int j = 0; j < S3_r[i].second.size() - 1; j++)
		{
			int node_id1 = S3_r[i].second[j],
				node_id2 = S3_r[i].second[j + 1],
				beta_id = 0;
			double hr = all_nodes[node_id2].r - all_nodes[node_id1].r,
				be = beta(all_nodes[node_id1].r, all_nodes[node_id1].z, beta_id);
			if (flag)
				GetM_Loc_dim2_r(all_nodes[node_id1].r, hr, M_loc);
			// mult M_loc on beta
			for (int k = 0; k < M_loc.size(); k++)
			{
				for (int l = 0; flag && l < M_loc[k].size(); l++)
					M_loc[k][l] *= be;
				b_loc[k] *= be;
			}
			if (flag)
				// add local to MS
				AddLocal_dim2(MS.ia, MS.ja, MS.di, MS.al, MS.au, M_loc, node_id1, node_id2);
			// get MS.b
			Getb_Loc_dim2_r(all_nodes[node_id1].r, all_nodes[node_id1].z, hr, b_loc, s3_id);
			MS.b.vect[node_id1] += b_loc[0];
			MS.b.vect[node_id2] += b_loc[1];
		}
	}
	NS2 = S3_z.size();
	for (int i = 0; i < NS2; i++)
	{
		int s3_id = S3_z[i].first;
		for (int j = 0; j < S3_z[i].second.size() - 1; j++)
		{
			int node_id1 = S3_z[i].second[j],
				node_id2 = S3_z[i].second[j + 1],
				beta_id = 0;
			double hz = all_nodes[node_id2].z - all_nodes[node_id1].z,
				be = beta(all_nodes[node_id1].r, all_nodes[node_id1].z, beta_id);
			if (flag)
				GetM_Loc_dim2_z(all_nodes[node_id1].z, hz, M_loc);
			// mult M_loc on beta
			for (int k = 0; k < M_loc.size(); k++)
			{
				for (int l = 0; flag && l < M_loc[k].size(); l++)
					M_loc[k][l] *= be;
				b_loc[k] *= be;
			}
			if (flag)
				// add local to MS
				AddLocal_dim2(MS.ia, MS.ja, MS.di, MS.al, MS.au, M_loc, node_id1, node_id2);
			// get MS.b
			Getb_Loc_dim2_z(all_nodes[node_id1].r, all_nodes[node_id1].z, hz, b_loc, s3_id);
			MS.b.vect[node_id1] += b_loc[0];
			MS.b.vect[node_id2] += b_loc[1];
		}
	}
	return 0;
}



int main()
{
	//LU_solver LU;
	//vector<int> 
	//	ia_t = {0, 0, 1, 2, 3, 5},
	//	ja_t = {0, 1, 1, 0, 2};
	//vector<double> 
	//	au_t = {-1, -2, -3, -4, -5},
	//	al_t = { 1, 2, 3, 4, 5 },
	//	di_t = {11, 22, 33, 44, 55},
	//	b_t = { -11, 2, 103, 182, 289 },
	//	x(5);

	//LU.Init(ia_t, ja_t, au_t, al_t, di_t, b_t);
	//LU.Calc(x); // right x = {1, 2, 3, 4, 5} if ja_t = {0, 1, 1, 0, 1}
	//b_t = { 47, 76, 107, 100, 95 };
	//LU.Recalc(b_t, x); // right x = {5, 4, 3, 2, 1} if ja_t = {0, 1, 1, 0, 1}

	ClearFolder("L:\\Мое\\курсач\\курсач\\output\\");
	//Node2D from(5, 0), to(7, 10);
	//func_2coord_to_1 test1 = test;
	gauss.Init();
	//double res = gauss.IntegrateFunc(from, to, test1);

	Make_grid2(""); // for creating tests
	Create_time_grid();

	Input();

	//std::vector<std::vector<double>> result(time_grid.size() - 1); // for tests

	MyMatrix M, G, A, MS;
	GeneratePortrait(M, all_nodes.size(), all_elems.size());
	G.ia = M.ia;
	G.ja = M.ja;
	M.au.resize(M.ja.size());
	M.al.resize(M.ja.size());
	M.N = all_nodes.size();
	M.di.resize(M.N);
	MS.ia = M.ia;
	MS.ja = M.ja;
	MS.au.resize(MS.ja.size());
	MS.al.resize(MS.ja.size());
	MS.N = all_nodes.size();
	MS.di.resize(MS.N);
	MS.b.Size(M.N);
	G.au.resize(G.ja.size());
	G.al.resize(G.ja.size());
	G.N = all_nodes.size();
	G.di.resize(G.N);


	std::vector<std::vector<double>> M_loc(4), G_loc(4), MG_loc(4);
	for (int i = 0; i < 4; i++)
	{
		M_loc[i].resize(4);
		G_loc[i].resize(4);
		MG_loc[i].resize(4);

	}
	std::vector<double>
		b_loc(4);

	// собираем матрицы M, G
	for (int i = 0; i < all_elems.size(); i++)
	{
		Get_Loc(M_loc, G_loc, MG_loc, i);
		AddLocal(M.ia, M.ja, M.di, M.al, M.au, M_loc, i);
		AddLocal(G.ia, G.ja, G.di, G.al, G.au, G_loc, i);
		AddLocal(G.ia, G.ja, G.di, G.al, G.au, MG_loc, i);
	}

	std::ofstream out;
	//out.imbue(std::locale("Russian"));
	out.precision(15);

	out.open("v.txt");
	vector<double> v(2);
	double r, z;
	int num = 0;
	out << all_elems.size() << '\n';
	for (int i = 0; i < all_elems.size(); i++)
	{
		r = (all_nodes[all_elems[i].node_loc[0]].r + all_nodes[all_elems[i].node_loc[3]].r) / 2;
		z = (all_nodes[all_elems[i].node_loc[0]].z + all_nodes[all_elems[i].node_loc[3]].z) / 2;
		num = GetVectorV(v,r,z,z_max_abs, r_max_abs);
		out << v[0] << " " << v[1] << " " << num << endl;
	}
	out.close();
	out.clear();


	double
		dt01 = 0,
		dt02 = 0,
		dt03 = 0,
		dt12 = 0,
		dt13 = 0,
		dt23 = 0;
	bool change_matrix;
	// явная
	/*for (i_t = 3; i_t < time_grid.size(); i_t++)
	{
		change_matrix = false;
		if (dt01 != time_grid[i_t] - time_grid[i_t - 1] ||
			dt02 != time_grid[i_t] - time_grid[i_t - 2] ||
			dt03 != time_grid[i_t] - time_grid[i_t - 3] ||
			dt12 != time_grid[i_t - 1] - time_grid[i_t - 2] ||
			dt13 != time_grid[i_t - 1] - time_grid[i_t - 3] ||
			dt23 != time_grid[i_t - 2] - time_grid[i_t - 3])
			change_matrix = true;

		dt01 = time_grid[i_t] - time_grid[i_t - 1];
		dt02 = time_grid[i_t] - time_grid[i_t - 2];
		dt03 = time_grid[i_t] - time_grid[i_t - 3];
		dt12 = time_grid[i_t - 1] - time_grid[i_t - 2];
		dt13 = time_grid[i_t - 1] - time_grid[i_t - 3];
		dt23 = time_grid[i_t - 2] - time_grid[i_t - 3];
		//пересобрать матрицу, если необходимо
		if (change_matrix)
		{
			A = M * ((dt01 * dt02 + dt01 * dt03 + dt02 * dt03) / (dt01 * dt02 * dt03));
			A.b.Size(G.N);
		}
		// пересобрать вектор правой части
		// обнуляет вектор
		for (int i = 0; i < A.b.vect.size(); i++)
		{
			A.b.vect[i] = 0;
			MS.b.vect[i] = 0;
		}
		i_t--; //по-идее, b с предыдущего времени, но работает только на текущем
		// собирает вектор
		for (int i = 0; i < all_elems.size(); i++)
		{
			Get_Loc_b(b_loc, i);
			AddLocal_b(A.b.vect, b_loc, i);
		}
		i_t++; //вернуть счетчик, если уменьшали
		MyVector temp;
		temp.Size(all_nodes.size());
		M.Ax(q1, temp);
		A.b = A.b + temp * ((dt01 * dt02) / (dt03 * dt13 * dt23));
		M.Ax(q2, temp);
		A.b = A.b + temp * ((-dt01 * dt03) / (dt02 * dt12 * dt23));
		M.Ax(q3, temp);
		A.b = A.b + temp * ((dt02 * dt03) / (dt01 * dt12 * dt13));
		G.Ax(q3, temp);
		A.b = A.b + temp * (-1);

		// учесть краевые
		//Set_S2(MS);
		//Set_S3(MS, change_matrix); //!!!!!!!!!!!!!!! матрица А не меняется
		//A = A + MS;
		//A.b = A.b + MS.b;
		SetS1(A.ia, A.ja, A.di, A.al, A.au, A.b.vect);

		//if (change_matrix) // если меняли матрицу, сбросить учтенные в А краевые
		//{
		//    std::fill(MS.al.begin(), MS.al.end(), 0);
		//    std::fill(MS.au.begin(), MS.au.end(), 0);
		//    std::fill(MS.di.begin(), MS.di.end(), 0);
		//}


		// решить СЛАУ
		Solver slau(A);
		slau.CGM_LU();
		slau.getx0(q4.vect);

		// вывести ответ на временном слое
		//out << "time = " << ";" << time_grid[i_t] << "\n";
		//for (int i = 0; i < all_nodes.size(); i++)
		//{
		//    out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q4.vect[i] << "\n";
		//}
		result[i_t - 3] = q4.vect; // for tests
		// сменить слой
		q1.vect.swap(q2.vect);
		q2.vect.swap(q3.vect);
		q3.vect.swap(q4.vect);
	}*/

	//Solver slau;
	//slau.InitMemory(all_nodes.size(), G.ja.size());
	// 
	LU_solver solver_LU;
	// неявная
	for (i_t = 3; i_t < time_grid.size(); i_t++)
	{
		cout << " i_t " << i_t << " time " << time_grid[i_t] << endl;
		change_matrix = false;
		if (dt01 != time_grid[i_t] - time_grid[i_t - 1] ||
			dt02 != time_grid[i_t] - time_grid[i_t - 2] ||
			dt03 != time_grid[i_t] - time_grid[i_t - 3] ||
			dt12 != time_grid[i_t - 1] - time_grid[i_t - 2] ||
			dt13 != time_grid[i_t - 1] - time_grid[i_t - 3] ||
			dt23 != time_grid[i_t - 2] - time_grid[i_t - 3])
			change_matrix = true;

		dt01 = time_grid[i_t] - time_grid[i_t - 1];
		dt02 = time_grid[i_t] - time_grid[i_t - 2];
		dt03 = time_grid[i_t] - time_grid[i_t - 3];
		dt12 = time_grid[i_t - 1] - time_grid[i_t - 2];
		dt13 = time_grid[i_t - 1] - time_grid[i_t - 3];
		dt23 = time_grid[i_t - 2] - time_grid[i_t - 3];

		double sum = 0;
		//пересобрать матрицу, если необходимо
		if (change_matrix)
		{
			A = G;
			A.b.Size(G.N);
			A = A + M * ((dt01 * dt02 + dt01 * dt03 + dt02 * dt03) / (dt01 * dt02 * dt03));

		}
		//sum += (dt01 * dt02 + dt01 * dt03 + dt02 * dt03) / (dt01 * dt02 * dt03);
		// пересобрать вектор правой части
		for (int i = 0; i < A.b.vect.size(); i++)
		{
			A.b.vect[i] = 0;
			MS.b.vect[i] = 0;
		}



		for (int i = 0; i < all_elems.size(); i++)
		{
			Get_Loc_b(b_loc, i);
			AddLocal_b(A.b.vect, b_loc, i);
		}
		MyVector temp;
		temp.Size(all_nodes.size());
		M.Ax(q1, temp);
		A.b = A.b + temp * ((dt01 * dt02) / (dt03 * dt13 * dt23));
		M.Ax(q2, temp);
		A.b = A.b + temp * ((-dt01 * dt03) / (dt02 * dt12 * dt23));
		M.Ax(q3, temp);
		A.b = A.b + temp * ((dt02 * dt03) / (dt01 * dt12 * dt13));

		// учесть краевые
		Set_S2(MS);
		Set_S3(MS, change_matrix); //!!!!!!!!!!!!!!! матрица А не меняется
		A = A + MS;
		A.b = A.b + MS.b;
		SetS1(A.ia, A.ja, A.di, A.al, A.au, A.b.vect);
		if (change_matrix) // если меняли матрицу, сбросить учтенные в А краевые
		{
			std::fill(MS.al.begin(), MS.al.end(), 0);
			std::fill(MS.au.begin(), MS.au.end(), 0);
			std::fill(MS.di.begin(), MS.di.end(), 0);
		}


		// решить СЛАУ
		//Solver slau(A);
		/*if (i_t == 3 || change_matrix)
		{
			slau.Clear(A);
		}
		else
		{
			slau.Clear(A.b.vect);
		}*/
		//slau.CGM_LU();
		//slau.getx0(q4.vect);

		if (i_t == 3 || change_matrix)
		{
			solver_LU.Init(A.ia, A.ja, A.au, A.al, A.di, A.b.vect);
			solver_LU.Calc(q4.vect);
		}
		else
		{
			solver_LU.Recalc(A.b.vect, q4.vect);
		}

		// вывести ответ на временном слое
		//out << "time = " << ";" << time_grid[i_t] << "\n";
		//for (int i = 0; i < all_nodes.size(); i++)
		//{
		//    out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q4.vect[i] << "\n";
		//}
		//result[i_t - 3] = q4.vect; // for tests
		out.open("output\\time_" + NumberToString(i_t) + ".txt");
		for (int i = 0; i < all_nodes.size(); i++)
		{
			out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q4.vect[i] << "\n"; // вывод в каждом узле
		}
		out.close();
		out.clear();
		out.open("output_elem\\time_" + NumberToString(i_t) + ".txt");
		// вывод в центре элемента
		for (int i = 0; i < all_elems.size(); i++)
		{
			node from = all_nodes[all_elems[i].node_loc[0]];
			node to = all_nodes[all_elems[i].node_loc[3]];
			node centre;
			centre.r = (from.r + to.r) / 2.0;
			centre.z = (from.z + to.z) / 2.0;

			double res = 0;
			for (int j = 0; j < 4; j++)
			{
				int node_id = all_elems[i].node_loc[j];
				res += q4.vect[node_id] * basis_1D[j % 2](from.r, to.r, centre.r) * basis_1D[j / 2](from.z, to.z, centre.z);
			}
			out << res << '\n';
		}
		out.close();
		out.clear();
		// сменить слой
		q1.vect.swap(q2.vect);
		q2.vect.swap(q3.vect);
		q3.vect.swap(q4.vect);
	}
	

	/*for (i_t = 1; i_t < time_grid.size(); i_t++)
	{
		cout << " i_t " << i_t << " time " << time_grid[i_t] << endl;
		change_matrix = false;
		if (dt01 != time_grid[i_t] - time_grid[i_t - 1])
			change_matrix = true;

		dt01 = time_grid[i_t] - time_grid[i_t - 1];
		double sum = 0;
		//пересобрать матрицу, если необходимо
		if (change_matrix)
		{
			A = G;
			A.b.Size(G.N);
			A = A + M * (1 / dt01);
		}

		// пересобрать вектор правой части
		for (int i = 0; i < A.b.vect.size(); i++)
		{
			A.b.vect[i] = 0;
			MS.b.vect[i] = 0;
		}
		for (int i = 0; i < all_elems.size(); i++)
		{
			Get_Loc_b(b_loc, i);
			AddLocal_b(A.b.vect, b_loc, i);
		}
		MyVector temp;
		temp.Size(all_nodes.size());
		M.Ax(q1, temp);
		A.b = A.b + temp * (1 / dt01);

		// учесть краевые
		Set_S2(MS);
		Set_S3(MS, change_matrix); //!!!!!!!!!!!!!!! матрица А не меняется
		A = A + MS;
		A.b = A.b + MS.b;
		SetS1(A.ia, A.ja, A.di, A.al, A.au, A.b.vect);
		if (change_matrix) // если меняли матрицу, сбросить учтенные в А краевые
		{
			std::fill(MS.al.begin(), MS.al.end(), 0);
			std::fill(MS.au.begin(), MS.au.end(), 0);
			std::fill(MS.di.begin(), MS.di.end(), 0);
		}



		// решить СЛАУ
		Solver slau(A);
		slau.CGM_LU();
		//slau.LOS_LU();
		slau.getx0(q2.vect);
		out.open("output\\time_" + NumberToString(i_t) + ".txt");
		for (int i = 0; i < all_nodes.size(); i++)
		{
			out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q2.vect[i] << "\n";
		}
		out.close();
		out.clear();
		// вывести ответ на временном слое
	   // out << "time = " << ";" << time_grid[i_t] << "\n";
		//for (int i = 0; i < all_nodes.size(); i++)
	   // {
	   //    out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q2.vect[i] << "\n";
	   // }
		result[i_t - 1] = q2.vect; // for tests
		//// сменить слой
		q1.vect.swap(q2.vect);

		//result[i_t - 3] = q4.vect; // for tests
	   // q1.vect.swap(q2.vect);
		//q2.vect.swap(q3.vect);
	   // q3.vect.swap(q4.vect);

	}
	*/
	// for tests
	/*bool outflag = false;
	for (int i = 0; i < all_nodes.size(); i++)
	{
		if ((int)(all_nodes[i].r * 10) % 10 == 0)
		{
			for (int j = 0; j < time_grid.size() - 3; j++)
			{
				if ((int)(all_nodes[i].z * 10) % 10 == 0)
				{
					out << result[j][i] << "\t";
					outflag = true;
				}
			}
			if (outflag)
			{
				out << "\n";
				outflag = false;
			}
		}
	}*/

	return 0;
}
