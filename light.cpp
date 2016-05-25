#include "light.h"

Color operator* (const Point3D &x, const Color &y)
{
    return Color(x.x*y.r, x.y*y.g, x.z*y.b);
}

Color operator+ (const Color &x, const Color &y)
{
    return Color(x.r+y.r, x.g+y.g, x.b+y.b);
}

Color operator* (const float &k, const Color &x)
{
    return Color(k*x.r, k*x.g, k*x.b);
}

float operator* (const Color &x, const Color &y)
{
    return x.r*y.r+x.g*y.g+x.b*y.b;
}
