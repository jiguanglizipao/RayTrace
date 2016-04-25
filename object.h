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
    float times;
    Point3D loc, Ks, Kt;
    bool readfile(std::string fi = "", float _times = 1.0, Point3D _loc = Point3D());
};

#endif
