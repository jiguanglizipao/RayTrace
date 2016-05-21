#ifndef POLYGON_H
#define POLYGON_H

#include "point.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include "matrix.h"

struct Polygon
{
    Point3D points3d[3], n;
    double d;
    Polygon(const Point3D _points3d[])
    {
        for(int i=0;i<3;i++)points3d[i] = _points3d[i];
        const Point3D &a = points3d[0], &b=points3d[1], &c=points3d[2];
        n = Point3D((b.y*c.z+c.y*a.z+a.y*b.z-b.y*a.z-a.y*c.z-c.y*b.z),
                    (a.x*c.z+b.x*a.z+c.x*b.z-c.x*a.z-b.x*c.z-a.x*b.z),
                    (a.x*b.y+b.x*c.y+c.x*a.y-c.x*b.y-b.x*a.y-a.x*c.y));
        d = -(a.x*b.y*c.z+b.x*c.y*a.z+c.x*a.y*b.z-c.x*b.y*a.z-b.x*a.y*c.z-a.x*c.y*b.z);

        double l = sqrt(n*n);
        n.norm(), d/=l;
    }

    double intersect(const Ray & r) const;

};

#endif
