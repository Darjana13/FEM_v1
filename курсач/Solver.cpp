#include "Solver.h"

Solver::Solver(int size)
{
    N = size;
    A.di.resize(N);
    A.ia.resize(N + 1);
    A.ja.resize((N * N - N) / 2);
    A.au.resize((N * N - N) / 2);
    A.al.resize((N * N - N) / 2);
    A.b.Size(N);
    A.N = N;

    A.ia[0] = 0;
    A.ia[1] = 0;
    A.di[0] = 1;
    A.b.vect[0] += 1;

    for (int i = 1; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            A.au[A.ia[i] + j] = 1. / (i + j + 1); // ò.ê. (i + 1) + (i - k + 1 + j) - 1, k = i
            A.b.vect[j] += (i + 1) * A.au[A.ia[i] + j];
            A.al[A.ia[i] + j] = 1. / (i + j + 1);
            A.b.vect[i] += (j + 1) * A.al[A.ia[i] + j];
            A.ja[A.ia[i] + j] = j;
        }
        A.ia[i + 1] = A.ia[i] + i;
        A.di[i] = 1. / (i + i + 1);
        A.b.vect[i] += (i + 1) * A.di[i];
    }
    maxIter = 10000;
    eps = 1E-14;

    std::ofstream fout;
    fout.precision(16);
    fout.open("kuslau.txt");
    fout << N << " " << maxIter << " " << eps;
    fout.close();

    fout.open("di.txt");
    for (int i = 0; i < N; i++)
        fout << A.di[i] << " ";
    fout.close();

    fout.open("ig.txt");
    for (int i = 0; i <= N; i++)
        fout << A.ia[i] << " ";
    fout.close();

    fout.open("jg.txt");
    for (int i = 0; i < A.ia[N]; i++)
        fout << A.ja[i] << " ";
    fout.close();

    fout.open("ggu.txt");
    for (int i = 0; i < A.ia[N]; i++)
        fout << A.au[i] << " ";
    fout.close();

    fout.open("ggl.txt");
    for (int i = 0; i < A.ia[N]; i++)
        fout << A.al[i] << " ";
    fout.close();

    fout.open("pr.txt");
    for (int i = 0; i < N; i++)
        fout << A.b.vect[i] << " ";
    fout.close();

    x0.Size(N);
    r.Size(N);
    z.Size(N);
    p.Size(N);
    Ar.Size(N);
    y.Size(N);
    L.resize(A.ia[N]);
    D.resize(N);
    U.resize(A.ia[N]);
    normB = A.b.Norm();
    iter = 0;
    normR = 0;
}

Solver::Solver(std::string filename)
{
    std::ifstream in("kuslau.txt");
    in >> N >> maxIter >> eps;
    A.ReadMatrix(N);
    x0.Size(N);
    r.Size(N);
    z.Size(N);
    p.Size(N);
    Ar.Size(N);
    y.Size(N);
    L.resize(A.ia[N]);
    D.resize(N);
    U.resize(A.ia[N]);
    normB = A.b.Norm();
    iter = 0;
    normR = 0;
}

Solver::Solver(MyMatrix _A)
{
    N = _A.N;
    maxIter = 10000;
    eps = 1E-15;
    A = _A;
    x0.Size(N);
    r.Size(N);
    z.Size(N);
    p.Size(N);
    Ar.Size(N);
    y.Size(N);
    L.resize(A.ia[N]);
    D.resize(N);
    U.resize(A.ia[N]);
    normB = A.b.Norm();
    iter = 0;
    normR = 0;
}

void Solver::output(std::string filename)
{
    std::ofstream out(filename);
    out.imbue(std::locale("Russian"));
    out.precision(15);
    for (int i = 0; i < N; i++)
        out << x0.vect[i] << std::endl;
}

void Solver::getx0(std::vector<double>& x)
{   
    for (int i = 0; i < N; i++)
        x[i] = x0.vect[i];
}

void Solver::CGM_LU()
{
    std::cout.precision(15);

    FactLU(L, U, D);

    double r_r = 0, Az_z = 0;
    double a = 0, B = 0;

    A.Ax(x0, r);          // r0 = A*x0
    A.b - r;              // r0 = B - A*x0
    Direct(L, D, r, r);   // r0 = L^(-1) * (B - A*x0)
    Reverse(L, D, r, r);  // r0 = L^(-T) * L^(-1) * (B - A*x0)
    A.ATx(r, y);          // y0 = A^(T) * L^(-T) * L^(-1) * (B - A*x0)
    Direct(U, r, y);      // r0 = U-t * A^(T) * L^(-T) * L^(-1) * (B - A*x0)

    z = r;                // z0 = r0
    r_r = r * r;
    normR = sqrt(r_r) / normB;
    for (iter = 1; iter < maxIter + 1 && normR >= eps; iter++)
    {
        Reverse(U, y, z);    // y = U^(-1) * z
        A.Ax(y, p);          // p = A * U^(-1) * z
        Direct(L, D, p, p);  // p = L-1 * A * U^(-1) * z
        Reverse(L, D, p, p); // p = L-t * p
        A.ATx(p, Ar);        // Ar = At * P
        Direct(U, Ar, Ar);   // Ar = U-t * Ar
        Az_z = Ar * z;       // (Ar,z)
        a = r_r / Az_z;

        // x(k) = x(k-1) + z(k-1)*a(k-1)
        // r(k) = r(k-1) - AT*A*z(k-1)*a(k-1)
        for (int i = 0; i < N; i++)
        {
            x0.vect[i] = x0.vect[i] + z.vect[i] * a;
            r.vect[i] = r.vect[i] - Ar.vect[i] * a;
        }

        // B(k) = (r(k), r(k)) / (r(k-1), r(k-1))
        B = 1.0 / r_r;
        r_r = r * r;
        B *= r_r;

        // z(k) = r(k) + B(k)*z(k-1)
        for (int i = 0; i < A.N; i++)
        {
            z.vect[i] = r.vect[i] + z.vect[i] * B;
        }
        normR = sqrt(r_r) / normB;
        //std::cout << iter << ". " << normR << std::endl;
    }
    // x0 = U^(-1) * x0
    Reverse(U, x0, x0);
}

void Solver::LOS_LU()
{
    std::cout.precision(15);

    FactLU(L, U, D);

    double p_p = 0, p_r = 0, r_r = 0, Ar_p = 0;
    double a = 0, B = 0, eps2 = 1e-10;

    A.Ax(x0, y);         // y = A * x0
    A.b - y;		         // y = B - A * x0
    Direct(L, D, r, y);  // r0 = L^(-1) * (B - A * x0)

    Reverse(U, z, r);    // z0 = U^(-1) * r0

    A.Ax(z, y);          // y = A * z0
    Direct(L, D, p, y);  // p0 = L^(-1) * (A * z0) 
    r_r = r * r;
    normR = sqrt(r_r) / normB;
    for (iter = 1; iter < maxIter + 1 && normR >= eps; iter++)
    {
        p_p = p * p;
        p_r = p * r;
        a = p_r / p_p;

        // x(k) = x(k-1) + a(k) * z(k-1)
        // r(k) = r(k-1) - a(k) * p(k-1)
        for (int i = 0; i < N; i++)
        {
            x0.vect[i] = x0.vect[i] + z.vect[i] * a;
            r.vect[i] = r.vect[i] - p.vect[i] * a;
        }

        Reverse(U, y, r);      // y = U^(-1) * r(k)
        A.Ax(y, Ar);           // Ar = A * U^(-1) * r(k)
        Direct(L, D, Ar, Ar);  // Ar = L^(-1) * A * U^(-1) * r(k)
        Ar_p = Ar * p;         // (Ar, p)
        B = -(Ar_p / p_p);

        // z(k) = U^(-1) * r(k) + B(k) * z(k-1) 
        // p(k) = L^(-1) * A * U^(-1) * r(k) + B(k) * p(k-1)
        for (int i = 0; i < N; i++)
        {
            z.vect[i] = y.vect[i] + z.vect[i] * B;
            p.vect[i] = Ar.vect[i] + p.vect[i] * B;
        }

        if (r_r - (r_r - a * a * p_p) < eps2)
            r_r = r * r;
        else
            r_r = r_r - a * a * p_p;
        normR = sqrt(r_r) / normB;
        std::cout << iter << ". " << normR << std::endl;
    }
}

void Solver::FactLU(std::vector<double>& L, std::vector<double>& U, std::vector<double>& D)
{
    L = A.al;
    U = A.au;
    D = A.di;
    double l, u, d;
    for (int k = 0; k < N; k++)
    {
        d = 0;
        int i0 = A.ia[k], i1 = A.ia[k + 1];
        int i = i0;
        for (; i0 < i1; i0++)
        {
            l = 0;
            u = 0;
            int j0 = i, j1 = i0;
            for (; j0 < j1; j0++)
            {
                int t0 = A.ia[A.ja[i0]], t1 = A.ia[A.ja[i0] + 1];
                for (; t0 < t1; t0++)
                {
                    if (A.ja[j0] == A.ja[t0])
                    {
                        l += L[j0] * U[t0];
                        u += L[t0] * U[j0];
                    }
                }
            }
            L[i0] -= l;
            U[i0] -= u;
            U[i0] /= D[A.ja[i0]];
            d += L[i0] * U[i0];
        }
        D[k] -= d;
    }
}

// L*y = B
void Solver::Direct(std::vector<double>& L, std::vector<double>& D, MyVector& y, MyVector& b)
{
    y = b;
    for (int i = 0; i < N; i++)
    {
        double sum = 0;
        int k0 = A.ia[i], k1 = A.ia[i + 1];
        int j;
        for (; k0 < k1; k0++)
        {
            j = A.ja[k0];
            sum += y.vect[j] * L[k0];
        }
        double buf = y.vect[i] - sum;
        y.vect[i] = buf / D[i];
    }
}

// U^(T)*y = B
void Solver::Direct(std::vector<double>& L, MyVector& y, MyVector& b)
{
    y = b;
    for (int i = 0; i < N; i++)
    {
        double sum = 0;
        int k0 = A.ia[i], k1 = A.ia[i + 1];
        int j;
        for (; k0 < k1; k0++)
        {
            j = A.ja[k0];
            sum += y.vect[j] * L[k0];
        }
        y.vect[i] -= sum;
    }
}

// U*x = y
void Solver::Reverse(std::vector<double>& U, MyVector& x, MyVector& y)
{
    x = y;
    for (int i = N - 1; i >= 0; i--)
    {
        int k0 = A.ia[i], k1 = A.ia[i + 1];
        int j;
        for (; k0 < k1; k0++)
        {
            j = A.ja[k0];
            x.vect[j] -= x.vect[i] * U[k0];
        }
    }
}

// L^(T)*x = y
void Solver::Reverse(std::vector<double>& U, std::vector<double>& D, MyVector& x, MyVector& y)
{
    x = y;
    for (int i = N - 1; i >= 0; i--)
    {
        int k0 = A.ia[i], k1 = A.ia[i + 1];
        int j;
        x.vect[i] /= D[i];
        for (; k0 < k1; k0++)
        {
            j = A.ja[k0];
            x.vect[j] -= x.vect[i] * U[k0];
        }
    }
}
