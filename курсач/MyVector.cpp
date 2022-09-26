

#include "MyVector.h"

MyVector::MyVector()
{
}

void MyVector::Size(int N)
{
    vect.resize(N);
}

// read from filename
void MyVector::ReadVector(std::string filename)
{
    if (vect.size() < 1)
        return;
    std::ifstream in(filename);
    for (int i = 0; i < vect.size(); i++)
    {
        in >> vect[i];
    }
    in.close();
}

// this = a; a = a
MyVector& MyVector::operator=(const MyVector& a)
{
    if (this != &a)
        this->vect = a.vect;
    return *this;
}

// this = a * this
MyVector& MyVector::operator*(const double a)
{
    for (int i = 0; i < this->vect.size(); i++)
        this->vect[i] *= a;
    return *this;
}

// (this, a)
double MyVector::operator* (const MyVector& a)
{
    double res = 0;
    if (this->vect.size() != a.vect.size())
        return res;

    for (int i = 0; i < this->vect.size(); i++)
        res += this->vect[i] * a.vect[i];
    return res;
}

// a = this - a;
MyVector& MyVector::operator-(MyVector& a)
{
    if (this->vect.size() != a.vect.size())
    {
        return *this;
    }
    else
    {
        for (int i = 0; i < this->vect.size(); i++)
        {
            a.vect[i] = this->vect[i] - a.vect[i];
        }
        return a;
    }
}
// this = this + a;
MyVector& MyVector::operator+(MyVector& a)
{
    if (this->vect.size() != a.vect.size())
    {
        return *this;
    }
    else
    {
        for (int i = 0; i < this->vect.size(); i++)
        {
            this->vect[i] = this->vect[i] + a.vect[i];
        }
        return *this;
    }
}

// || this || 
double MyVector::Norm()
{
    return sqrt((*this) * (*this));
}
