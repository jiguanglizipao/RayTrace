#include "light.h"

Color operator* (const Point3D &x, const Color &y)
{
    return Color(x.x*y.r, x.y*y.g, x.z*y.b);
}
