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
