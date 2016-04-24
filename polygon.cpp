#include "point.h"
#include "polygon.h"

Polygon::State Polygon::checkInside(Point2D pos)
{
    int sum = 0;
    for(std::size_t i=0;i<points2d.size();i++){
        Point2D x = points2d[i], y = points2d[(i+1)%points2d.size()];
        if(onLine(x, y, pos))
            return online;
        else
        {
            x = x-pos, y=y-pos;
            int kx = x.getQuadrant(), ky = y.getQuadrant();
            if((4+ky-kx)%4 == 1)sum += 1;
            if((4+ky-kx)%4 == 2)sum += ((x^y)>0)?2:-2;
            if((4+ky-kx)%4 == 3)sum -= 1;
        }
    }
    if(sum%8 == 0)return outside;
    if(abs(sum)%8 == 4)return inside;
    return error;
}
