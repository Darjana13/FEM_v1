#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include "Solver.h"
#include "Generate.h"

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

std::vector<node> all_nodes;         // ??? ???? ? ??????? ?????????? ?????????
std::vector<element> all_elems;      // ??? ???????? ? ?????? ?????????? ?????????
std::vector<material> all_materials; // ??? ????????? ?? ????????
std::vector<std::pair<int,std::vector<int>>> S1;      // S1[i][j] ?? j-?? ???? ?????? ??????? 1 ????
std::vector<std::pair<int,std::vector<int>>> S2_r;    // ??????? ??????????? r ?? j-?? ???? ?????? ??????? 2 ????
std::vector<std::pair<int,std::vector<int>>> S2_z;    // ??????? ??????????? z ?? j-?? ???? ?????? ??????? 2 ????
std::vector<std::pair<int,std::vector<int>>> S3_r;    // ??????? ??????????? r ?? j-?? ???? ?????? ??????? 3 ????
std::vector<std::pair<int,std::vector<int>>> S3_z;    // ??????? ??????????? z ?? j-?? ???? ?????? ??????? 3 ????

std::vector<double> time_grid; // ????? ???????
int i_t = 0; // ??????? ????????? ????

double gamma(double r, double z, int gam_id) // ???????? ????? ?? ??????? gam_id 
{
    switch (gam_id)
    {
    case 0:
        return 1;
    default:
        std::cout << "can't find gamma ? " << gam_id << "\n";
        break;
    }
}
double beta(double r, double z, int beta_id) // ???????, ??? ???? ????? ??????????
{
    switch (beta_id)
    {
    case 0:
        return 1;
    default:
        std::cout << "can't find gamma ? " << beta_id << "\n";
        break;
    }
}

double func_f(double r, double z, int f_id) // ???????? f ?? ??????? f_id 
{
    double t = time_grid[i_t];
    switch (f_id)
    {
    case 0:
        return   -4 - 5

            ;
    default:
        std::cout << "can't find f ? " << f_id << "\n";
        break;
    }
}

double func_S(double r, double z, int s_id) // ???????? ???????? S ?? ??????? f_id
{
    double t = time_grid[i_t];
    switch (s_id)
    {
    case 0:
        return r * r + 3 * z - 5 * t

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
        std::cout << "can't find S ? " << s_id << "\n";
        break;
    }
}

int Input() // ?????? ??????
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

    // ??? ????????????
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
    // ????? ??? ????????????

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
    std::vector<std::vector<double>>& G_loc) // ????????? ????????? G
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
    std::vector<std::vector<double>>& M_loc) // ??????????? ????????? ?
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
    std::vector<double>& b_loc, int f_id) // ????????? ?????????? b
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
    int el_id) // ????????? ????????? ??????? ?
{
    element el = all_elems[el_id];
    double hr = all_nodes[el.node_loc[1]].r - all_nodes[el.node_loc[0]].r,
        hz = all_nodes[el.node_loc[2]].z - all_nodes[el.node_loc[0]].z;

    // ???????? ?
    GetM_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, all_materials[el.mater].gamma_id, hr, hz, M_loc);
    
    // ???????? G
    GetG_Loc(all_nodes[el.node_loc[0]].r, all_materials[el.mater].lambda, hr, hz, G_loc); // A_loc = G_loc     
    return 0;
}
int Get_Loc_b(std::vector<double>& b_loc,
    int el_id) // ???????? ?????? ?????? ?????
{
    element el = all_elems[el_id];
    double hr = all_nodes[el.node_loc[1]].r - all_nodes[el.node_loc[0]].r,
        hz = all_nodes[el.node_loc[2]].z - all_nodes[el.node_loc[0]].z;
    // ???????? f
    Getb_Loc(all_nodes[el.node_loc[0]].r, all_nodes[el.node_loc[0]].z, hr, hz, b_loc, el.f_id);
    return 0;
}
int GeneratePortrait(MyMatrix &A, 
    int N, int Kel) // ????????? ????????
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
    // ???????? ????????? A, b  ? ?????????? ????
{ 
    std::vector<int> L = all_elems[el_id].node_loc;
    int n_loc = all_elems[el_id].node_loc.size(); // ??????????? ????????? ???????
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
    // ???????? ????????? b  ? ?????????? ????
{
    std::vector<int> L = all_elems[el_id].node_loc;
    int k = all_elems[el_id].node_loc.size(); // ??????????? ????????? ???????
    for (int i = 0; i < k; i++)
    {
        b[L[i]] += b_loc[i];
    }
    return 0;
}

int SetS1(std::vector<int>& ia, std::vector<int>& ja, std::vector<double>& di,
    std::vector<double>& al, std::vector<double>& au,
    std::vector<double>& b) // ???? ?????? ???????
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
    std::vector<std::vector<double>>& M_loc) // ????????? ??????? ??? S3_r
{
    
    M_loc[0][0] = hr / 6 * (2 * rp + hr / 2 );
    M_loc[0][1] = hr / 6 * (rp + hr / 2);
    M_loc[1][0] = hr / 6 * (rp + hr / 2);
    M_loc[1][1] = hr / 6 * (2 * rp + 3 * hr / 2);

    return 0;
}
int Getb_Loc_dim2_r(double rp, double zs, double hr,
    std::vector<double>& b_loc, int f_id) // ????????? ?????????? b ??? S ???????????? r
{
    double f1 = func_S(rp, zs, f_id),
        f2 = func_S(rp + hr, zs, f_id);
    b_loc[0] = f1 * (hr / 6 * (2 * rp + hr / 2)) + f2 * (hr / 6 * (rp + hr / 2));
    b_loc[1] = f1 * (hr / 6 * (rp + hr / 2)) + f2 * (hr / 6 * (2 * rp + 3 * hr / 2));
    return 0;
}
double GetM_Loc_dim2_z(double zp, double hz,
    std::vector<std::vector<double>>& M_loc) // ????????? ??????? ??? S3_z
{
    M_loc[0][0] = hz / 3;
    M_loc[0][1] = hz / 6;
    M_loc[1][0] = hz / 6;
    M_loc[1][1] = hz / 3;

    return 0;
}
int Getb_Loc_dim2_z(double rp, double zs, double hz,
    std::vector<double>& b_loc, int s_id) // ????????? ?????????? b ??? S ???????????? z
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
    int node1, int node2) // ???????? ????????? ??????? ? ?????????? ??? ??????????? ??????
{
    std::vector<int> L(2);
    L[0] = node1;
    L[1] = node2;
    int n_loc = 2; // ??????????? ????????? ???????
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
int Set_S2(MyMatrix &MS) // ???? ?????? ??????? MS.b - ?????? ?????? ?? ???????
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
int Set_S3(MyMatrix& MS, bool flag) // MS.b - ?????? ?????? ?? ???????
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

    Input();

    std::vector<std::vector<double>> result(time_grid.size()-3); // for tests

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
    

    std::vector<std::vector<double>> M_loc(4), G_loc(4);
    for (int i = 0; i < 4; i++)
    {
        M_loc[i].resize(4);
        G_loc[i].resize(4);
    }
    std::vector<double>
        b_loc(4);

    // ???????? ??????? M, G
    for (int i = 0; i < all_elems.size(); i++)
    {
        Get_Loc(M_loc, G_loc, i);
        AddLocal(M.ia, M.ja, M.di, M.al, M.au, M_loc, i);
        AddLocal(G.ia, G.ja, G.di, G.al, G.au, G_loc, i);
    }

    std::ofstream out("result.txt");
    //out.imbue(std::locale("Russian"));
    out.precision(15);

    double 
        dt01 = 0,
        dt02 = 0,
        dt03 = 0,
        dt12 = 0,
        dt13 = 0,
        dt23 = 0;
    bool change_matrix;
    // ?????
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
        //??????????? ???????, ???? ??????????
        if (change_matrix)
        {
            A = M * ((dt01 * dt02 + dt01 * dt03 + dt02 * dt03) / (dt01 * dt02 * dt03));
            A.b.Size(G.N);
        }
        // ??????????? ?????? ?????? ?????
        // ???????? ??????
        for (int i = 0; i < A.b.vect.size(); i++)
        {
            A.b.vect[i] = 0;
            MS.b.vect[i] = 0;
        }
        i_t--; //??-????, b ? ??????????? ???????, ?? ???????? ?????? ?? ???????
        // ???????? ??????
        for (int i = 0; i < all_elems.size(); i++)
        {
            Get_Loc_b(b_loc, i);
            AddLocal_b(A.b.vect, b_loc, i);
        }
        i_t++; //??????? ???????, ???? ?????????
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

        // ?????? ???????
        //Set_S2(MS);
        //Set_S3(MS, change_matrix); //!!!!!!!!!!!!!!! ??????? ? ?? ????????
        //A = A + MS;
        //A.b = A.b + MS.b;
        SetS1(A.ia, A.ja, A.di, A.al, A.au, A.b.vect);
        
        //if (change_matrix) // ???? ?????? ???????, ???????? ???????? ? ? ???????
        //{
        //    std::fill(MS.al.begin(), MS.al.end(), 0);
        //    std::fill(MS.au.begin(), MS.au.end(), 0);
        //    std::fill(MS.di.begin(), MS.di.end(), 0);
        //}


        // ?????? ????
        Solver slau(A);
        slau.CGM_LU();
        slau.getx0(q4.vect);

        // ??????? ????? ?? ????????? ????
        //out << "time = " << ";" << time_grid[i_t] << "\n";
        //for (int i = 0; i < all_nodes.size(); i++)
        //{
        //    out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q4.vect[i] << "\n";
        //}
        result[i_t - 3] = q4.vect; // for tests
        // ??????? ????
        q1.vect.swap(q2.vect);
        q2.vect.swap(q3.vect);
        q3.vect.swap(q4.vect);
    }*/

    // ???????
    for (i_t = 3; i_t < time_grid.size(); i_t++)
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

        double sum = 0;
        //??????????? ???????, ???? ??????????
        if (change_matrix)
        {         
            A = G;
            A.b.Size(G.N);
            A = A + M * ((dt01 * dt02 + dt01 * dt03 + dt02 * dt03) / (dt01 * dt02 * dt03));  
          
        }
        //sum += (dt01 * dt02 + dt01 * dt03 + dt02 * dt03) / (dt01 * dt02 * dt03);
        // ??????????? ?????? ?????? ?????
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

        //sum += (dt01 * dt02) / (dt03 * dt13 * dt23);
        //sum += (-dt01 * dt03) / (dt02 * dt12 * dt23);
        //sum += (dt02 * dt03) / (dt01 * dt12 * dt13);
        //std::cout << sum << "\n";

        // ?????? ???????
        Set_S2(MS);
        Set_S3(MS, change_matrix); //!!!!!!!!!!!!!!! ??????? ? ?? ????????
        A = A + MS;
        A.b = A.b + MS.b;
        SetS1(A.ia, A.ja, A.di, A.al, A.au, A.b.vect);
        if (change_matrix) // ???? ?????? ???????, ???????? ???????? ? ? ???????
        {
            std::fill(MS.al.begin(), MS.al.end(), 0);
            std::fill(MS.au.begin(), MS.au.end(), 0);
            std::fill(MS.di.begin(), MS.di.end(), 0);
        }


        // ?????? ????
        Solver slau(A);
        slau.CGM_LU();
        slau.getx0(q4.vect);

        // ??????? ????? ?? ????????? ????
        //out << "time = " << ";" << time_grid[i_t] << "\n";
        //for (int i = 0; i < all_nodes.size(); i++)
        //{
        //    out << all_nodes[i].r << "\t" << all_nodes[i].z << "\t" << q4.vect[i] << "\n";
        //}
        result[i_t - 3] = q4.vect; // for tests
        // ??????? ????
        q1.vect.swap(q2.vect);
        q2.vect.swap(q3.vect);
        q3.vect.swap(q4.vect);
    }

    // for tests
    bool outflag = false;
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
    }

    return 0;
}
