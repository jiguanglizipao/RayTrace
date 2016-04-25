#ifndef LIGHT_H
#define LIGHT_H

#include "point.h"

struct Color
{
    float r, g, b;
    Color(float _r = 0, float _g = 0, float _b = 0)
        :r(_r), g(_g), b(_b)
    {
    }
    Color operator+ (const Color &x)
    {
        return Color(r+x.r, g+x.g, b+x.b);
    }

    Color operator* (const float &k)
    {
        return Color(k*r, k*g, k*b);
    }

    float operator* (const Color &x)
    {
        return x.r*r+x.g*g+x.b*b;
    }
};

Color operator* (const Point3D &x, const Color &y);

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
