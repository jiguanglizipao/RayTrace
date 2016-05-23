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
    std::vector<cv::Mat*> tex;
    std::vector<std::string> name;
    double times;
    Point3D loc, col, lig;
    RType type;
    bool readfile(std::string fi = "", double _times = 1.0, Point3D _loc = Point3D(), Point3D rotate = Point3D());
};

#endif
