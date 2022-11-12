#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <set>
#include "Solver.h"
#include "Generate.h"
using namespace std;

typedef // создаем новый прототип (в данном случае указатель на функцию)
double // возвращаемое значение (такое же как в функциях)
(*basis_func) // имя прототипа (в коде употребляется без звездочки)
(double, double, double); // список параметров (такое же как в функциях)

double linear1(double x1, double x2, double x) { return (x2 - x) / (x2 - x1); }
double linear2(double x1, double x2, double x) { return (x - x1) / (x2 - x1); }
double dlinear1(double x1, double x2, double x) { return (-1.0) / (x2 - x1); }
double dlinear2(double x1, double x2, double x) { return (1.0) / (x2 - x1); }


MyVector q1, q2, q3, q4;

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
std::vector<std::pair<int,std::vector<int>>> S1;      // S1[i][j] на j-ом узле заданы краевые 1 рода
std::vector<std::pair<int,std::vector<int>>> S2_r;    // граница параллельна r на j-ом узле заданы краевые 2 рода
std::vector<std::pair<int,std::vector<int>>> S2_z;    // граница параллельна z на j-ом узле заданы краевые 2 рода
std::vector<std::pair<int,std::vector<int>>> S3_r;    // граница параллельна r на j-ом узле заданы краевые 3 рода
std::vector<std::pair<int,std::vector<int>>> S3_z;    // граница параллельна z на j-ом узле заданы краевые 3 рода

std::vector<double> time_grid; // сетка времени
int i_t = 0; // текущий временной слой

double gamma(double r, double z, int gam_id) // значение гамма по индексу gam_id 
{
    switch (gam_id)
    {
    case 0:
        // !!! сюда написать коэффициент gamma
        return 1;
    default:
        std::cout << "can't find gamma № " << gam_id << "\n";
        break;
    }
    return 1;

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
    return 1;

}

double func_f(double r, double z, int f_id) // значение f по индексу f_id 
{
    double t = 0;
    if(time_grid.size() != 0)
        t = time_grid[i_t];
    switch (f_id)
    {
    case 0:
        return  //5
            // !!! сюда написать -div (lambda grad u*) + gamma u*, (посмотреть, как раскрывается div grad в цилиндрических координатах)
            // u* указать в func_S 
            // lambda - первое число в файле material.txt
            // gamma в функции gamma
            11.0 * cos(3.0 * r + z) + 3.0 / r * sin(3.0 * r + z)

            ;
    default:
        std::cout << "can't find f № " << f_id << "\n";
        break;
    }
}

double func_S(double r, double z, int s_id) // значение краевого S по индексу f_id
{
    double t = 0;
    if (time_grid.size() != 0)
        t = time_grid[i_t];
    switch (s_id)
    {
    case 0:
        return //5
            // !!! сюда написать истинную функцию u*
            cos(3.0 * r + z)

            ;
    case 1:// 2_z
        return 2 *r*r
            ;
    case 2: // 2_r
        return 3 
            ;
    case 3: // 3_z
        return 1 / beta(r,z,0) * 2 * r * r + r*r +3*z - 5*t
            ;
    case 4: // 3_r
        return 1 / beta(r, z, 0) * 3 + r * r + 3 * z - 5 * t
            ;
    default:
        std::cout << "can't find S № " << s_id << "\n";
        break;
    }
}

int Input_no_time() // чтение данных
{
    int N, Nmat, Kel, NS1, Ntime, NS;
    std::ifstream in;

    in.open("info.txt");
    //in >> N >> Nmat >> Kel >> NS1;
    Nmat = 1;
    NS1 = 1;
    in.close();

    in.open("rz.txt");
    in >> N;
    all_nodes.resize(N);
    for (int i = 0; i < N; i++)
    {
        in >> all_nodes[i].r >> all_nodes[i].z;
    }
    in.close();

    /*in.open("time.txt");
    in >> Ntime;
    time_grid.resize(Ntime);
    for (int i = 0; i < Ntime; i++)
    {
        in >> time_grid[i];
    }
    in.close();

    // для тестирования
    std::ofstream temp("q0 q1 q2.txt");
    temp.precision(15);
    for (int k = 0; k < 3; k++)
    {
        i_t = k;
        for (int i = 0; i < N; i++)
        {
            temp << func_S(all_nodes[i].r, all_nodes[i].z, 0) << " ";
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
    }*/
    q4.Size(N);
    //in.close();

    in.open("S1.txt");
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


int Input() // чтение данных
{
    int N, Nmat, Kel, NS1, Ntime, NS;
    std::ifstream in;

    in.open("info.txt");
    in >> N >> Nmat >> Kel >> NS1;
    in.close();

    in.open("rz.txt");
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
    std::ofstream temp("q0 q1 q2.txt");
    temp.precision(15);
    for (int k = 0; k < 3; k++)
    {
        i_t = k;
        for (int i = 0; i < N; i++)
        {
            temp << func_S(all_nodes[i].r, all_nodes[i].z, 0) << " ";
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
    all_materials.resize(Nmat);
    for (int i = 0; i < Nmat; i++)
    {
        in >> all_materials[i].lambda >> all_materials[i].gamma_id;
    }
    in.close();

    in.open("elem.txt");
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

int Get_Loc(std::vector<std::vector<double>>& M_loc, std::vector<std::vector<double>>& G_loc,
    int el_id) // получение локальной матрицы А
{
    element el = all_elems[el_id];
    double hr = all_nodes[el.node_loc[1]].r - all_nodes[el.node_loc[0]].r,
        hz = all_nodes[el.node_loc[2]].z - all_nodes[el.node_loc[0]].z;

    // получить М
    GetM_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, all_materials[el.mater].gamma_id, hr, hz, M_loc);
    
    // получить G
    GetG_Loc(all_nodes[el.node_loc[0]].r, all_materials[el.mater].lambda, hr, hz, G_loc); // A_loc = G_loc     
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
int GeneratePortrait(MyMatrix &A, 
    int N, int Kel) // генерация портрета
{
    int N_loc = 4;
    A.ia.resize(N + 1);
   
    vector<set<int>> map(N);
    int k = 0;
    for (auto elem : all_elems)
    {
        for (auto i : elem.node_loc)
            for (auto j : elem.node_loc)
            {
                if (i > j)
                  map[i].insert(j);
            }
    }
    A.ia[0] = 0;
    for (int i = 0; i < N; i++)
    {
        A.ia[i + 1] = A.ia[i] + map[i].size();
    }

    A.ja.resize(A.ia[N]);
    for (int i = 0; i < N; i++)
    {
        set <int> ::iterator it = map[i].begin();
        for (int j = 0; it != map[i].end(); it++, j++)
        {
            A.ja[A.ia[i] + j] = *(it);
        }
    }

    return 0;

    /*std::vector<int>* ia = &A.ia,
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
    return 0;*/
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
    
    M_loc[0][0] = hr / 6 * (2 * rp + hr / 2 );
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
    double f1 = func_S(rp, zs     , s_id);
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
int Set_S2(MyMatrix &MS) // учет вторых краевых MS.b - вектор вклада от краевых
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
            if(flag)
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
    Make_grid(""); // for creating tests
    //Create_time_grid();

    //Input();
    Input_no_time();

    //std::vector<std::vector<double>> result(time_grid.size()-3); // for tests

    MyMatrix M, G, A, MS;
    GeneratePortrait(M, all_nodes.size(), all_elems.size());
    G.ia = M.ia;
    G.ja = M.ja;
    A.ia = M.ia;
    A.ja = M.ja;
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
    A.au.resize(G.ja.size());
    A.al.resize(G.ja.size());
    A.N = all_nodes.size();
    A.di.resize(G.N);
    A.b.vect.resize(G.N);

    std::vector<std::vector<double>> M_loc(4), G_loc(4);
    for (int i = 0; i < 4; i++)
    {
        M_loc[i].resize(4);
        G_loc[i].resize(4);
    }
    std::vector<double>
        b_loc(4);

    // собираем матрицы M, G
    for (int i = 0; i < all_elems.size(); i++)
    {
        Get_Loc(M_loc, G_loc, i);
        AddLocal(A.ia, A.ja, A.di, A.al, A.au, M_loc, i);
        AddLocal(A.ia, A.ja, A.di, A.al, A.au, G_loc, i);
        Get_Loc_b(b_loc, i);
        AddLocal_b(A.b.vect, b_loc, i);
    }
    //A = A + G;
    Set_S2(MS);
    Set_S3(MS, true); //!!!!!!!!!!!!!!! матрица А не меняется
    A = A + MS;
    A.b = A.b + MS.b;
    SetS1(A.ia, A.ja, A.di, A.al, A.au, A.b.vect);
    std::fill(MS.al.begin(), MS.al.end(), 0);
    std::fill(MS.au.begin(), MS.au.end(), 0);
    std::fill(MS.di.begin(), MS.di.end(), 0);
    
    Solver slau(A);
    slau.CGM_LU();
    slau.getx0(q4.vect);


    // !!! тут вывод
    ofstream res("Point.txt"); // !!! этот файл надо будет перетащить в папку со сплайном
    ofstream file("f(r,z).txt"); // !!! в этом файле пишется r, z, f, df/dr, df/dz, ddf/(drdz), надо сравнить с аналогичной штукой в сплайне

    vector<basis_func> basis1D = { linear1, linear2 };
    vector<basis_func> dbasis1D = { dlinear1, dlinear2 };

    int p_count = all_elems.size();
    double tmp, r, z, f, df_dx, df_dy, df_dx_dy;
    double r1, r2, z1, z2;
    res << p_count << endl;
    double sum_residual = 0;
    for (int el_i = 0; el_i < all_elems.size(); el_i++)
    {
        r1 = all_nodes[all_elems[el_i].node_loc[0]].r;
        r2 = all_nodes[all_elems[el_i].node_loc[1]].r;
        z1 = all_nodes[all_elems[el_i].node_loc[0]].z;
        z2 = all_nodes[all_elems[el_i].node_loc[2]].z;

        r = (r1 + r2) / 2;
        z = (z1 + z2) / 2;

        f = 0;
        df_dx = 0;
        df_dy = 0;
        df_dx_dy = 0;
        for (int i = 0; i < 4; i++)
        {
            f += q4.vect[all_elems[el_i].node_loc[i]] * basis1D[i%2](r1, r2, r) * basis1D[i/2](z1, z2, z);
            df_dx += q4.vect[all_elems[el_i].node_loc[i]] * dbasis1D[i % 2](r1, r2, r) * basis1D[i / 2](z1, z2, z);
            df_dy += q4.vect[all_elems[el_i].node_loc[i]] * basis1D[i % 2](r1, r2, r) * dbasis1D[i / 2](z1, z2, z);
            df_dx_dy += q4.vect[all_elems[el_i].node_loc[i]] * dbasis1D[i % 2](r1, r2, r) * dbasis1D[i / 2](z1, z2, z);
        }
        res << r << '\t' << z << '\t' << f << '\t' << 1 << endl;
        file << r << '\t' << z << '\t' << f << '\t' << df_dx << '\t' << df_dy << '\t' << df_dx_dy << endl;
        //check << r << '\t' << z << '\t' << f << '\t' << func_S(r, z, 0) << '\t' << abs(f - func_S(r, z, 0)) / abs(func_S(r, z, 0)) << endl;
        sum_residual += abs(f - func_S(r, z, 0)) / abs(func_S(r, z, 0));
    }
    cout << sum_residual / p_count << endl; // !!! в консоль выведется средняя относительная погрешность решения
    res.close();
    res.clear();
    file.close();
    file.clear();



    return 0;
}
