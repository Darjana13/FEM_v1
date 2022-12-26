#include "LU_solver.h"

void LU_solver::Init(vector<int> &_ia, vector<int>& _ja, vector<double>& _au, vector<double>& _al, vector<double>& _di, vector<double>& _b)
{
    N = _di.size();
    di.resize(N);
    for (int i = 0; i < N; i++)
        di[i] = _di[i];

    al.clear();
    au.clear();
    al.reserve(_ia[N] * 2);
    au.reserve(_ia[N] * 2);
    ia.resize(N + 1);

    int col = 0, count = 0;
    for (int i = 0; i < N; i++)
    {
        ia[i] = count;
        col = _ja[_ia[i]];
        for (int j = _ia[i]; j < _ia[i + 1]; j++)
        {
            // i - row, _ja[j] - col
            if (col == _ja[j])
            {
                al.push_back(_al[j]);
                au.push_back(_au[j]);
                count++;
                col++;
            }
            else
            {
                while (col != _ja[j])
                {
                    al.push_back(0.0);
                    au.push_back(0.0);
                    count++;
                    col++;
                }
                al.push_back(_al[j]);
                au.push_back(_au[j]);
                count++;
                col++;
            }
        }
        while (col != i)
        {
            al.push_back(0.0);
            au.push_back(0.0);
            count++;
            col++;
        }
    }
    ia[N] = count;
    b = _b;
    y.resize(N);

    return;
}


void LU_solver::Read(string path)
{
    ifstream fin;
    fin.open("IN.txt");

    fin >> N;
    di.resize(N);
    for (int i = 0; i < N; i++)
        fin >> di[i];

    ia.resize(N + 1);
    for (int i = 0; i <= N; i++)
        fin >> ia[i];

    au.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
        fin >> au[i];

    al.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
        fin >> al[i];

    b.resize(N);
    for (int i = 0; i < N; i++)
        fin >> b[i];

    fin.close();

    return;
}

void LU_solver::forward_stroke() // Прямой ход Ly = F
{
    double res = 0;
    y = b;
    for (int i = 0; i < N; i++)
    {
        int count = i - (ia[i + 1] - ia[i]);
        res = 0;
        for (int j = ia[i]; j < ia[i + 1]; j++, count++)
        {
            res += y[count] * al[j];
        }
        y[i] = (y[i] - res) / di[i];
    }
}

void LU_solver::revers()
{
    for (int i = N - 1; i >= 0; i--)
    {
        int j = i - (ia[i + 1] - ia[i]);
        for (int k = ia[i]; k < ia[i + 1]; k++, j++)
        {
            y[j] -= au[k] * y[i];
        }
    }

}

bool LU_solver::LU()
{

    for (int i = 1; i < N; i++)
    {
        int j0 = i - (ia[i + 1] - ia[i]);
        for (int ii = ia[i]; ii < ia[i + 1]; ii++)
        {
            int j = ii - ia[i] + j0;
            double sum_l = 0, sum_u = 0;
            if (ia[j] < ia[j + 1])
            {
                int j0j = j - (ia[j + 1] - ia[j]);
                int jjbeg = j0 < j0j ? j0j : j0; // max (j0, j0j)
                int jjend = j < i - 1 ? j : i - 1; // min (j, i - 1)
                for (int k = 0; k < jjend - jjbeg; k++)
                {
                    int ind_prev = ia[j] + jjbeg - j0j + k;
                    int ind_now = ia[i] + jjbeg - j0 + k;
                    sum_l += au[ind_prev] * al[ind_now];
                    sum_u += au[ind_now] * al[ind_prev];
                }
            }
            al[ii] -= sum_l;
            au[ii] -= sum_u;
            if (abs(di[j]) < eps) // matrix hasn't LU
            {
                // cout << "di[" << j << "] = " << di[j] << endl;
                return false;
            }
            au[ii] /= di[j];
            di[i] -= al[ii] * au[ii];
        }
    }
    au_LU = au;
    al_LU = al;
    di_LU = di;

    return true;
}

void LU_solver::Recalc(vector<double>& _b, vector<double>& _x)
{
    b = _b;

    au = au_LU;
    al = al_LU;
    di = di_LU;

    forward_stroke();
    revers();

    _x = y;

    return;
}

void LU_solver::Calc(vector<double>& _x)
{
    LU();

    forward_stroke();
    revers();

    _x = y;

    return;
}
