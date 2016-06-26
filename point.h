#ifndef POINT_H
#define POINT_H
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>

const double eps = 1e-5;

enum RType { DIFF, SPEC, REFR, BDR };	// material types, used in radiance()

struct Point3D
{
    double x, y, z;

    __device__ __host__ Point3D(double t=0)
    {
        x = t;
        y = t;
        z = t;
    }

    __device__ __host__ Point3D(double x_, double y_, double z_)
    {
        x = x_;
        y = y_;
        z = z_;
    }

    __device__ __host__ Point3D operator+(const Point3D & b) const
    {
        return Point3D(x + b.x, y + b.y, z + b.z);
    }

    __device__ __host__ Point3D operator-(const Point3D & b) const
    {
        return Point3D(x - b.x, y - b.y, z - b.z);
    }

    __device__ __host__ Point3D operator*(double b) const
    {
        return Point3D(x * b, y * b, z * b);
    }

    __device__ __host__ Point3D mult(const Point3D & b) const
    {
        return Point3D(x * b.x, y * b.y, z * b.z);
    }

    __device__ __host__ Point3D & norm() {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    __device__ __host__ double operator*(const Point3D & b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }

    __device__ __host__ Point3D operator%(const Point3D & b) const{
        return Point3D(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

struct Ray
{
    Point3D o, d;
    __device__ __host__ Ray(Point3D o_, Point3D d_)
        :o(o_), d(d_)
    {
    }
};


#endif
