#ifndef POINT_H
#define POINT_H

const float eps = 1e-7;
struct Point3D
{
    float x, y, z;
    Point3D(float _x = 0.0, float _y = 0.0, float _z = 0.0)
        :x(_x), y(_y), z(_z)
    {
    }
};

struct Point2D
{
    float x, y;
    Point2D(float _x = 0.0, float _y = 0.0)
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
Point2D operator*(const float &k, const Point2D &x);
float operator*(const Point2D &x, const Point2D &y);
float operator^(const Point2D &x, const Point2D &y);
float operator|(const Point2D &x, const Point2D &y);

Point3D operator+(const Point3D &x, const Point3D &y);
Point3D operator-(const Point3D &x, const Point3D &y);
float operator*(const Point3D &x, const Point3D &y);
Point3D operator*(const float &k, const Point3D &x);
Point3D operator^(const Point3D &x, const Point3D &y);


#endif
