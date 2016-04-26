#ifndef POLYGON_H
#define POLYGON_H

#include "point.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include "matrix.h"

struct Polygon
{
    std::vector<Point2D> points2d;
    std::vector<Point3D> points3d;
    Matrix rotate, rotate_r;
    Point3D n;
    double d;
    enum State
    {
        outside, inside
    };

    Polygon(const std::vector<Point3D> &_points3d)
        :points3d(_points3d)
    {
        const Point3D &a = points3d[0], &b=points3d[1], &c=points3d[2];
        n = Point3D((b.y*c.z+c.y*a.z+a.y*b.z-b.y*a.z-a.y*c.z-c.y*b.z),
                    (a.x*c.z+b.x*a.z+c.x*b.z-c.x*a.z-b.x*c.z-a.x*b.z),
                    (a.x*b.y+b.x*c.y+c.x*a.y-c.x*b.y-b.x*a.y-a.x*c.y));
        d = -(a.x*b.y*c.z+b.x*c.y*a.z+c.x*a.y*b.z-c.x*b.y*a.z-b.x*a.y*c.z-a.x*c.y*b.z);

        double l = sqrt(n*n);
        n=(1/l)*n, d/=l;

        rotate = Matrix::Ry(sqrt(n.y*n.y+n.z*n.z)/sqrt(n*n), -n.x/sqrt(n*n))*Matrix::Rx(n.z/sqrt(n.y*n.y+n.z*n.z), n.y/sqrt(n.z*n.z+n.y*n.y))*Matrix::T(-a.x, -a.y, -a.z);
        rotate_r = Matrix::T(a.x, a.y, a.z)*Matrix::Rx(n.z/sqrt(n.y*n.y+n.z*n.z), -n.y/sqrt(n.z*n.z+n.y*n.y))*Matrix::Ry(sqrt(n.y*n.y+n.z*n.z)/sqrt(n*n), n.x/sqrt(n*n));

        points2d.clear();
        for(std::size_t i=0;i<points3d.size();i++)
        {
            Point3D t = rotate*points3d[i];
            points2d.push_back(Point2D(t.x, t.y));
        }
    }

    bool SameSide(Point3D a, Point3D b, Point3D c, Point3D p) const;

    State checkInside(Point3D pos) const;

};

#endif
