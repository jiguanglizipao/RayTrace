#ifndef LIGHT_H
#define LIGHT_H

#include "point.h"

struct Color
{
    double r, g, b;
    Color(double _r = 0, double _g = 0, double _b = 0)
        :r(_r), g(_g), b(_b)
    {
    }
};

Color operator* (const Point3D &x, const Color &y);
Color operator+ (const Color &x, const Color &y);
Color operator* (const double &k, const Color &x);
double operator* (const Color &x, const Color &y);
struct Light
{
    Point3D loc;
    Color col;
    Light(Point3D _loc = Point3D(), Color _col = Color())
        :loc(_loc), col(_col)
    {
    }
};

#endif // LIGHT_H
