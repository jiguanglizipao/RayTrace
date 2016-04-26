#include "point.h"
#include "polygon.h"
#include <cstdio>

bool Polygon::SameSide(Point3D a, Point3D b, Point3D c, Point3D p) const
{
    return (((b-a)^(c-a))*((b-a)^(p-a))) > -eps;
}

Polygon::State Polygon::checkInside(Point3D p) const
{
    Point3D a=points3d[0], b=points3d[1], c=points3d[2];
    if(SameSide(a, b, c, p) && SameSide(b, c, a, p) && SameSide(c, a, b, p))
        return inside;
    return outside;
}
