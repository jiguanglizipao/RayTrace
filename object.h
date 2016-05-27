#ifndef OBJECT_H
#define OBJECT_H
#include "point.h"
#include "polygon.h"
#include <vector>
#include <string>

struct Object
{
    std::vector<Polygon> polys;
    std::vector<Point3D> points, cols;
    std::vector<unsigned char*> tex;
    std::vector<int> sizex, sizey;
    std::vector<std::string> name;
    float times;
    Point3D loc, col, lig;
    RType type;
    bool readfile(std::string fi = "", float _times = 1.0, Point3D _loc = Point3D(), Point3D rotate = Point3D());
};

#endif
