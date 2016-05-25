#include "sphere.h"

bool sphere_intersect(Sphere spheres[], int n, const Ray & r, int &id, double &t)
{
    double d;
    t = 1e20;
    for (int i = int (n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < 1e10;
}
