#ifndef SPHERE_H
#define SPHERE_H
#include "point.h"

struct Sphere
{
    float rad;
    Point3D pos, lig, col;
    RType type;
    __device__ __host__ Sphere(float rad_, Point3D pos_, Point3D lig_, Point3D col_, RType type_)
        :rad(rad_), pos(pos_), lig(lig_), col(col_), type(type_) {
    }
    __device__ __host__ float intersect(const Ray & r) const
    {
        Point3D op = pos - r.o;
        float t, b = op * r.d, det = b * b - op * op + rad * rad;
        if (det < 0)return 1e20;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 1e20);
    }
};

bool sphere_intersect(const Sphere spheres[], int n, const Ray & r, int &id, float &t);

#endif
