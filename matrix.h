#ifndef MATRIX_H
#define MATRIX_H
#include <cstring>
#include <cmath>
#include "point.h"

struct Matrix
{
    double a[4][4], n, m;
    Matrix(int _n, int _m)
        :n(_n), m(_m)
    {
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)a[i][j]=0;
    }

    static Matrix T(double x, double y, double z)
    {
        Matrix ans(4, 4);
        ans.a[0][0]=ans.a[1][1]=ans.a[2][2]=ans.a[3][3] = 1;
        ans.a[0][3]=x, ans.a[1][3]=y, ans.a[2][3]=z;
        return ans;
    }

    static Matrix Rx(double cost, double sint)
    {
        Matrix ans(4, 4);
        ans.a[0][0]=ans.a[3][3]=1;
        ans.a[1][1]=ans.a[2][2]=cost;
        ans.a[1][2]=ans.a[2][1]=sint;
        ans.a[1][2]*=-1;
        return ans;
    }

    static Matrix Ry(double cost, double sint)
    {
        Matrix ans(4, 4);
        ans.a[1][1]=ans.a[3][3]=1;
        ans.a[0][0]=ans.a[2][2]=cost;
        ans.a[0][2]=ans.a[2][0]=sint;
        ans.a[2][0]*=-1;
        return ans;
    }

    static Matrix Rz(double cost, double sint)
    {
        Matrix ans(4, 4);
        ans.a[2][2]=ans.a[3][3]=1;
        ans.a[0][0]=ans.a[1][1]=cost;
        ans.a[0][1]=ans.a[1][0]=sint;
        ans.a[0][1]*=-1;
        return ans;
    }
    static Matrix RotateMatrix(Point3D s, Point3D v, double theta);
};

Matrix operator+ (Matrix a, Matrix b);
Matrix operator- (Matrix a, Matrix b);
Matrix operator* (Matrix a, Matrix b);
Point3D operator* (Matrix b, Point3D a);

#endif // MATRIX_H
