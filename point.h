#ifndef POINT_H
#define POINT_H

const float eps = 1e-10;
struct Point3D                                                          //3维点对
{
    float x, y, z;
    Point3D(double _x = 0.0, double _y = 0.0, double _z = 0.0)
        :x(_x), y(_y), z(_z)
    {
    }
};

struct Point2D                                                          //2维点对          
{
    float x, y;
    Point2D(float _x = 0.0, float _y = 0.0)
        :x(_x), y(_y)
    {
    }
    int getQuadrant()                                                   //获得点的象限
    {
        if(x>-eps && y>-eps)return 1;
        if(x<-eps && y>-eps)return 2;
        if(x<-eps && y<-eps)return 3;
        if(x>-eps && y<-eps)return 4;
    }
};

Point2D operator-(const Point2D &x, const Point2D &y);                  //2维点（向量）相减
Point2D operator+(const Point2D &x, const Point2D &y);                  //2维点（向量）相加
Point2D operator*(const float &k, const Point2D &x);                   //向量数乘
double operator*(const Point2D &x, const Point2D &y);                   //向量点乘
double operator^(const Point2D &x, const Point2D &y);                   //向量叉乘模
double operator|(const Point2D &x, const Point2D &y);                   //向量长度

#endif
