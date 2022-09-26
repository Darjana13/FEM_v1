#include "MyMatrix.h"

MyMatrix::MyMatrix(void)
{
}

void MyMatrix::ReadMatrix(int size)
{
    N = size;
    std::ifstream in;

    in.open("ig.txt");
    ia.resize(N + 1);
    for (int i = 0; i < N + 1; i++)
    {
        in >> ia[i];
    }
    in.close();
    if (ia[0])
        for (int i = 0; i < N + 1; i++)
        {
            ia[i]--;
        }

    in.open("jg.txt");
    ja.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
    {
        in >> ja[i];
    }
    in.close();
    if (ja[0])
        for (int i = 0; i < ia[N]; i++)
        {
            ja[i]--;
        }

    in.open("di.txt");
    di.resize(N);
    for (int i = 0; i < N; i++)
    {
        in >> di[i];
    }
    in.close();

    in.open("ggu.txt");
    au.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
    {
        in >> au[i];
    }
    in.close();

    in.open("ggl.txt");
    al.resize(ia[N]);
    for (int i = 0; i < ia[N]; i++)
    {
        in >> al[i];
    }
    in.close();

    b.Size(N);
    b.ReadVector("pr.txt");
}

// y = Ax
void MyMatrix::Ax(MyVector& x, MyVector& y)
{
    for (int i = 0; i < N; i++)
    {
        y.vect[i] = di[i] * x.vect[i];
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {
            int k = ja[j];
            y.vect[i] += al[j] * x.vect[k];
            y.vect[k] += au[j] * x.vect[i];
        }
    }
}

// y = Ax
void MyMatrix::Ax(std::vector<double>& x, std::vector<double>& y)
{
    for (int i = 0; i < N; i++)
    {
        y[i] = di[i] * x[i];
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {
            int k = ja[j];
            y[i] += al[j] * x[k];
            y[k] += au[j] * x[i];
        }
    }
}

// y = A^(T)x
void MyMatrix::ATx(MyVector& x, MyVector& y)
{
    for (int i = 0; i < N; i++)
    {
        y.vect[i] = di[i] * x.vect[i];
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {
            int k = ja[j];
            y.vect[i] += au[j] * x.vect[k];
            y.vect[k] += al[j] * x.vect[i];
        }
    }
}

MyMatrix& MyMatrix::operator+ (MyMatrix B)
{
    if (N != B.N)
    {
        std::cout << "A и B разного размера\n";
        return *this;
    }
    for (int i = 0; i < N; i++)
    {
        this->di[i] += B.di[i];
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {
            int k = ja[j];
            if (k != B.ja[j])
            {
                std::cout << "A и B имеют разные портреты\n";
                return *this;
            }
            this->al[j] += B.al[j];
            this->au[j] += B.au[j];
        }
    }
    return *this;
}

MyMatrix MyMatrix::operator* (const double a)
{
    MyMatrix C = *this;
    for (int i = 0; i < N; i++)
    {
        C.di[i] *= a;
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {
            C.al[j] *= a;
            C.au[j] *= a;
        }
    }
    return C;
}

MyMatrix& MyMatrix::operator=(const MyMatrix& B)
{
    if (this != &B)
    {
        this->al = B.al;
        this->au = B.au;
        this->b = B.b;
        this->di = B.di;
        this->ia = B.ia;
        this->ja = B.ja;
        this->N = B.N;
    }
    return *this;
}

void MyMatrix::Show()
{
    std::vector<std::vector<double>> Matr(N);
    for (int i = 0; i < N; i++)
    {
        Matr[i].resize(N);
    }
    for (int i = 0; i < N; i++)
    {       
        for (int j = ia[i]; j < ia[i + 1]; j++)
        {

        }
    }
}