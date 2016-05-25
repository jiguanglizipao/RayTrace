#include "point.h"
#include "polygon.h"
#include <cstdio>

static inline unsigned char GetImageData(unsigned char *image, int sizex, int sizey, double x, double y, int c)
{
//    int t = (int(x*sizex+0.5)%sizex)*3*sizey+3*(int(y*sizey+0.5)%sizey)+c;
//    printf("%d %d\n", t, 3*sizex*sizey);
    return image[(int(x*sizex+0.5)%sizex)*3*sizey+3*(int(y*sizey+0.5)%sizey)+c];
}


double Polygon::intersect(const Ray & r, double &u, double &v) const
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

Point3D Polygon::getn(Point3D x) const
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

Point3D Polygon::getcol(double u, double v) const
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
