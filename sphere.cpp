#include "sphere.h"

bool sphere_intersect(const Sphere spheres[], int n, const Ray & r, int &id, float &t)
{
    float d;
    t = 1e20;
    for (int i = int (n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < 1e10;
}
