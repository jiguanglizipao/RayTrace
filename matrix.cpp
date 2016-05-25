#include "matrix.h"

Matrix operator+ (Matrix a, Matrix b)
{
    Matrix ans(a.n, a.m);
    for(int i=0;i<a.n;i++)
        for(int j=0;j<a.m;j++)
            ans.a[i][j]=a.a[i][j]+b.a[i][j];
    return ans;
}

Matrix operator- (Matrix a, Matrix b)
{
    Matrix ans(a.n, a.m);
    for(int i=0;i<a.n;i++)
        for(int j=0;j<a.m;j++)
            ans.a[i][j]=a.a[i][j]-b.a[i][j];
    return ans;
}

Matrix operator* (Matrix a, Matrix b)
{
    Matrix ans(a.n, b.m);
    for(int i=0;i<a.n;i++)
        for(int k=0;k<a.m;k++)
        for(int j=0;j<b.m;j++)
            ans.a[i][j]+=a.a[i][k]*b.a[k][j];
    return ans;
}

Point3D operator* (Matrix b, Point3D a)
{
    Matrix tmp(4, 1);
    tmp.a[0][0]=a.x, tmp.a[1][0]=a.y, tmp.a[2][0]=a.z, tmp.a[3][0]=1;
    Matrix ans = b*tmp;
    return Point3D(ans.a[0][0], ans.a[1][0], ans.a[2][0]);
}

Matrix Matrix::RotateMatrix(Point3D s, Point3D p, float theta)
{
    float u=p.x, v=p.y, w=p.z;
    Matrix a(4, 4);
    a.a[0][0]=u*u+(v*v+w*w)*cos(theta);
    a.a[0][1]=u*v*(1-cos(theta))-w*sin(theta);
    a.a[0][2]=u*w*(1-cos(theta))+v*sin(theta);
    a.a[0][3]=(s.x*(v*v+w*w)-u*(s.y*v+s.z*w))*(1-cos(theta))+(s.y*w-s.z*v)*sin(theta);
    a.a[1][0]=u*v*(1-cos(theta))+w*sin(theta);
    a.a[1][1]=v*v+(u*u+w*w)*cos(theta);
    a.a[1][2]=v*w*(1-cos(theta))-u*sin(theta);
    a.a[1][3]=(s.y*(u*u+w*w)-v*(s.x*u+s.z*w))*(1-cos(theta))+(s.z*u-s.x*w)*sin(theta);
    a.a[2][0]=u*w*(1-cos(theta))-v*sin(theta);
    a.a[2][1]=v*w*(1-cos(theta))+u*sin(theta);
    a.a[2][2]=w*w+(u*u+v*v)*cos(theta);
    a.a[1][3]=(s.z*(u*u+v*v)-w*(s.x*u+s.y*v))*(1-cos(theta))+(s.x*v-s.y*u)*sin(theta);
    a.a[3][0]=a.a[3][1]=a.a[3][2]=0;a.a[3][3]=1;
    return a;
}
