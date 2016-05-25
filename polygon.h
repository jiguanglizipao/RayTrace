#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include "point.h"
#include "matrix.h"

struct Polygon
{
    Point3D points3d[3], n[3], ts[3], lig, col, nc;
    unsigned char *tex;
    int sizex, sizey;
    RType type;
    __host__ Polygon(const Point3D points3d_[], const Point3D n_[], const Point3D ts_[], Point3D lig_, Point3D col_, RType type_, unsigned char *tex_, int sizex_, int sizey_)
        :lig(lig_), col(col_), type(type_), tex(tex_), sizex(sizex_), sizey(sizey_)
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

    __device__ __host__ unsigned char GetImageData(unsigned char *image, int sizex, int sizey, double x, double y, int c) const
    {
    //    int t = (int(x*sizex+0.5)%sizex)*3*sizey+3*(int(y*sizey+0.5)%sizey)+c;
    //    printf("%d %d\n", t, 3*sizex*sizey);
        return image[(int(x*sizex+0.5)%sizex)*3*sizey+3*(int(y*sizey+0.5)%sizey)+c];
    }

    __device__ __host__ double intersect(const Ray & r, double &u, double &v) const
    {
        Point3D E1 = points3d[1] - points3d[0];
        Point3D E2 = points3d[2] - points3d[0];
        Point3D P = r.d%E2;
        double det = E1*P, t;
        Point3D T;
        if(det > 0)
            T = r.o - points3d[0];
        else
            T = points3d[0] - r.o, det = -det;
        if(det < eps)return 1e20;
        u = T*P;
        if( u < 0 || u > det )return 1e20;
        Point3D Q = T%E1;
        v = r.d*Q;
        if(v < 0 || u + v > det )return 1e20;
        t = E2*Q;
        double fInvDet = 1.0 / det;
        t *= fInvDet;
        u *= fInvDet;
        v *= fInvDet;
        return t;
    }
    __device__ __host__ Point3D getn(Point3D x) const
    {
        Point3D v[3], t[3], ans;
        double s[3], sum=0;
        for(int i=0;i<3;i++)v[i] = (points3d[i]-x).norm();
        for(int i=0;i<3;i++)t[i] = v[i]%v[(i+1)%3];
        for(int i=0;i<3;i++)s[i] = sqrt(t[i]*t[i]),sum+=s[i];
        sum*=0.1;
        if(s[0]>sum || s[1]>sum || s[2]>sum)return nc;
        for(int i=0;i<3;i++)ans = ans + n[i]*s[(i+1)%3];
        return ans*double(1.0/(s[0]+s[1]+s[2]));
    }
    __device__ __host__ Point3D getcol(double u, double v) const
    {
        if(!tex)return col;
        Point3D s = (ts[1]-ts[0])*u+(ts[2]-ts[0])*v+ts[0];
        Point3D ans;
    //    if(s.x > 2 || s.y > 2)printf("%lf %lf\n", s.x, s.y);
        ans.x = GetImageData(tex, sizex, sizey, s.x, s.y, 2);
        ans.y = GetImageData(tex, sizex, sizey, s.x, s.y, 1);
        ans.z = GetImageData(tex, sizex, sizey, s.x, s.y, 0);
    //    printf("%lf %lf %lf\n", ans.x, ans.y, ans.z);

    //    ans.x = tex->at<cv::Vec3b>(s.x,s.y)[2];
    //    ans.y = tex->at<cv::Vec3b>(s.x,s.y)[1];
    //    ans.z = tex->at<cv::Vec3b>(s.x,s.y)[0];
        ans.x = pow(ans.x/256, 2.2);ans.y = pow(ans.y/256, 2.2);ans.z = pow(ans.z/256, 2.2);
        return ans;
    }

};

#endif
