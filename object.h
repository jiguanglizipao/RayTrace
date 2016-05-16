#ifndef OBJECT_H
#define OBJECT_H
#include "point.h"
#include "polygon.h"
#include <vector>
#include <string>

struct Object
{
    std::vector<Polygon> polys;
    std::vector<Point3D> points;
    double times;
    Point3D loc, Ks, Kt, Kds, Kdt, Ks1, Kt1, Ks2, Kt2;
    bool readfile(std::string fi = "", double _times = 1.0, Point3D _loc = Point3D(), Point3D rotate = Point3D());
};

#endif
