#ifndef POLYGON_H
#define POLYGON_H

#include "point.h"
#include <vector>
#include <algorithm>
#include <cmath>
struct Polygon
{
    std::vector<Point2D> points2d;                                                          //3维点信息
    std::vector<Point3D> points3d;                                                          //在x-y平面上投影的2维点信息
    Point3D n;                                                                              //平面法向量
    float d;                                                                                //单位法向量组成的方程中的d
    enum State                                                                              //判断点是否在面上的状态类
    {
        online, outside, inside, error
    };

    Polygon(Point3D _n, float _d, const std::vector<Point3D> &_points3d)
        :n(_n), d(_d), points3d(_points3d)
    {
        points2d.clear();
        for(std::size_t i=0;i<points3d.size();i++)
            points2d.push_back(Point2D(points3d[i].x, points3d[i].y));
    }

    Polygon(const std::vector<Point3D> &_points3d)                                          //在没有法向量等时计算法向量
        :points3d(_points3d)
    {
        points2d.clear();
        for(std::size_t i=0;i<points3d.size();i++)
            points2d.push_back(Point2D(points3d[i].x, points3d[i].y));
        const Point3D &a = points3d[0], &b=points3d[1], &c=points3d[2];
        n = Point3D((b.y*c.z+c.y*a.z+a.y*b.z-b.y*a.z-a.y*c.z-c.y*b.z),
                    (a.x*c.z+b.x*a.z+c.x*b.z-c.x*a.z-b.x*c.z-a.x*b.z),
                    (a.x*b.y+b.x*c.y+c.x*a.y-c.x*b.y-b.x*a.y-a.x*c.y));
        d = -(a.x*b.y*c.z+b.x*c.y*a.z+c.x*a.y*b.z-c.x*b.y*a.z-b.x*a.y*c.z-a.x*c.y*b.z);
    }

    double disptoseg(Point2D a, Point2D b, Point2D p)                                       //点到线段距离
    {
        Point2D t(p.x+a.y-b.y, p.y+b.x-a.x);
        if((((a-p)^(t-p))*((b-p)^(t-p))) > eps)
            return ((p|a)<(p|b))?(p|a):(p|b);
        return fabs((p-b)^(a-b))/(a|b);
    }

    bool onLine(Point2D a, Point2D b, Point2D pos)                                          //判断点是否在线上
    {
        return disptoseg(a, b, pos)<0.5;
    }

    State checkInside(Point2D pos);                                                         //判断点是否在面内

    double getDepth(Point2D pos)                                                            //获得一个坐标在该面上的深度
    {
        if(n.z < eps)return -1e30;
        return -(d+n.x*pos.x+n.y*pos.y)/n.z;
    }
};

#endif
