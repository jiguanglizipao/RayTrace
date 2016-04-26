#include "point.h"
#include <cmath>

Point2D operator-(const Point2D &x, const Point2D &y){
    return Point2D(x.x-y.x, x.y-y.y);
}

Point2D operator+(const Point2D &x, const Point2D &y){
    return Point2D(x.x+y.x, x.y+y.y);
}

Point2D operator*(const double &k, const Point2D &x){
    return Point2D(k*x.x, k*x.y);
}

double operator*(const Point2D &x, const Point2D &y){
    return x.x*y.x+x.y*y.y;
}

double operator^(const Point2D &x, const Point2D &y){
    return x.x*y.y-y.x*x.y;
}

double operator|(const Point2D &x, const Point2D &y){
    return sqrt((x.x-y.x)*(x.x-y.x)+(y.y-x.y)*(y.y-x.y));
}

Point3D operator+(const Point3D &x, const Point3D &y){
    return Point3D(x.x+y.x, x.y+y.y, x.z+y.z);
}

Point3D operator-(const Point3D &x, const Point3D &y){
    return Point3D(x.x-y.x, x.y-y.y, x.z-y.z);
}

double operator*(const Point3D &x, const Point3D &y){
    return x.x*y.x+x.y*y.y+x.z*y.z;
}

Point3D operator*(const double &k, const Point3D &x){
    return Point3D(k*x.x, k*x.y, k*x.z);
}

Point3D operator^(const Point3D &x, const Point3D &y){
    return Point3D(x.y*y.z-x.z*y.y, x.z*y.x-x.x*y.z, x.x*y.y-x.y*y.x);
}
