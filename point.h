#ifndef POINT_H
#define POINT_H

const double eps = 1e-5;
struct Point3D
{
    double x, y, z;
    Point3D(double _x = 0.0, double _y = 0.0, double _z = 0.0)
        :x(_x), y(_y), z(_z)
    {
    }
};

struct Point2D
{
    double x, y;
    Point2D(double _x = 0.0, double _y = 0.0)
        :x(_x), y(_y)
    {
    }
    int getQuadrant()
    {
        if(x>-eps && y>-eps)return 1;
        if(x<-eps && y>-eps)return 2;
        if(x<-eps && y<-eps)return 3;
        if(x>-eps && y<-eps)return 4;
    }
};

Point2D operator-(const Point2D &x, const Point2D &y);
Point2D operator+(const Point2D &x, const Point2D &y);
Point2D operator*(const double &k, const Point2D &x);
double operator*(const Point2D &x, const Point2D &y);
double operator^(const Point2D &x, const Point2D &y);
double operator|(const Point2D &x, const Point2D &y);

Point3D operator+(const Point3D &x, const Point3D &y);
Point3D operator-(const Point3D &x, const Point3D &y);
double operator*(const Point3D &x, const Point3D &y);
Point3D operator*(const double &k, const Point3D &x);
Point3D operator^(const Point3D &x, const Point3D &y);


#endif
