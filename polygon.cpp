#include "point.h"
#include "polygon.h"
#include <cstdio>

double Polygon::intersect(const Ray & r) const
{
    Point3D E1 = points3d[1] - points3d[0];
    Point3D E2 = points3d[2] - points3d[0];
    Point3D P = r.d%E2;
    double det = E1*P, u, v, t;
    Point3D T;
    if(det > 0)
        T = r.o - points3d[0];
    else
        T = points3d[0] - r.o, det = -det;
    if(det < eps)return 1e20;
    u = T*P;
    if( u < 0 || u > det )return 1e20;
    Point3D Q = T%E1;
    v = r.d*Q;
    if(v < 0 || u + v > det )return 1e20;
    t = E2*Q;
    double fInvDet = 1.0 / det;
    t *= fInvDet;
    u *= fInvDet;
    v *= fInvDet;
    return t;
}
