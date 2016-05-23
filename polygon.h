#ifndef POLYGON_H
#define POLYGON_H

#include "point.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cv.hpp>
//#if CV_VERSION_MAJOR == 3
//#include <opencv2/highgui.hpp>
//#else
//#include <opencv/highgui.h>
//#endif
#include "matrix.h"

struct Polygon
{
    Point3D points3d[3], n[3], ts[3], lig, col, nc;
    const cv::Mat *tex;
    Polygon(const Point3D points3d_[], const Point3D n_[], const Point3D ts_[], Point3D lig_, Point3D col_, const cv::Mat *tex_)
        :lig(lig_), col(col_), tex(tex_)
    {
        for(int i=0;i<3;i++)points3d[i] = points3d_[i];
        for(int i=0;i<3;i++)n[i] = n_[i]-points3d_[i];
        for(int i=0;i<3;i++)ts[i] = ts_[i];
        const Point3D &a = points3d[0], &b=points3d[1], &c=points3d[2];
        nc = Point3D((b.y*c.z+c.y*a.z+a.y*b.z-b.y*a.z-a.y*c.z-c.y*b.z),
                    (a.x*c.z+b.x*a.z+c.x*b.z-c.x*a.z-b.x*c.z-a.x*b.z),
                    (a.x*b.y+b.x*c.y+c.x*a.y-c.x*b.y-b.x*a.y-a.x*c.y));
        nc.norm();
    }

    double intersect(const Ray & r, double &u, double &v) const;
    Point3D getn(Point3D x) const;
    Point3D getcol(double u, double v) const;

};

#endif
