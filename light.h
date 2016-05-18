#ifndef LIGHT_H
#define LIGHT_H

#include "point.h"
#include "kdtree.h"

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
    double ma[3], mi[3];
    Color col;
    double check_aabb(Point3D view, Point3D ray)
    {
        double l = sqrt(ray*ray), s[3]={view.x, view.y, view.z}, v[3]={ray.x/l, ray.y/l, ray.z/l};
        Point3D pos;
        if (mi[0]<s[0]+eps && s[0]<ma[0]+eps && mi[1]<s[1]+eps && s[1]<ma[1]+eps && mi[2]<s[2]+eps && s[2]<ma[2]+eps)return 0;
        double dis = 1e15;
        for(int i=0;i<3;i++)
        {
            if (fabs(v[i]) > eps)
            {
                double t = ((v[i]>0?mi[i]:ma[i]) - s[i]) / v[i];
                if (t > eps)
                {
                    pos = view + t * ray;
                    double p[3]={pos.x, pos.y, pos.z};
                    if (mi[(i+1)%3]<p[(i+1)%3]+eps && p[(i+1)%3]<ma[(i+1)%3]+eps && mi[(i+2)%3]<p[(i+2)%3]+eps && p[(i+2)%3]<ma[(i+2)%3]+eps)dis = std::min(dis, fabs(t));
                }
            }
        }
        return dis;
    }
    Light(Color _col = Color())
        :col(_col)
    {
    }
};

#endif // LIGHT_H
