#include <cv.h>
#if CV_VERSION_MAJOR == 3
    #include <opencv2/highgui.hpp>
#else
    #include <opencv/highgui.h>
#endif
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cmath>
#include "point.h"
#include "polygon.h"
#include "object.h"
#include "kdtree.h"
#include "light.h"

using namespace std;
using namespace cv;

vector<Object> objs;
vector<Light> lights;
KdTree *kdtree;

void drawPixel(cv::Mat &image, int x, int y, Color color) {
    Color ans = color;
    if(ans.b > 1 || ans.r > 1 || ans.g > 1)
        printf("%lf %lf %lf\n", ans.r, ans.g, ans.b);
    ans.r = ans.r>1?1:ans.r;
    ans.g = ans.g>1?1:ans.g;
    ans.b = ans.b>1?1:ans.b;
    image.at<cv::Vec3b>(x,y)[2] = int(ans.r*255);
    image.at<cv::Vec3b>(x,y)[1] = int(ans.g*255);
    image.at<cv::Vec3b>(x,y)[0] = int(ans.b*255);
}

bool CheckAll(Point3D view, Point3D ray, int &no, int &nv, Point3D &p, double &dis)
{
    dis=1.1e10;
    double raydis = sqrt(ray*ray);
    for(size_t i=0;i<objs.size();i++)
    {
        for(size_t j=0;j<objs[i].polys.size();j++)
        {
            double k = -(objs[i].polys[j].d+objs[i].polys[j].n*view)/(objs[i].polys[j].n*ray);
            if(k < eps)continue;
            Point3D tmp = view+k*ray;
            if(objs[i].polys[j].checkInside(tmp) == Polygon::inside)
                if(k*raydis < dis)dis=k*raydis, no=i, nv=j, p=tmp;
        }
    }
    return dis < 1e10;
}

Color Ilocal(Point3D view, Point3D pos, int no, int nv)
{
    Color ans(0.0, 0.0, 0.0);
    for(size_t i=0;i<lights.size();i++)
    {
        int no1, nv1;
        Point3D p, L=lights[i].loc-pos, V=view-pos, tmp(0.0, 0.0, 0.0), n;
        if(L*objs[no].polys[nv].n < 0)n=-1*objs[no].polys[nv].n;else n=objs[no].polys[nv].n;
        double dis;
        kdtree->check(objs, pos, L, no1, nv1, p, dis);
        //CheckAll(pos, L, no1, nv1, p, dis);
        if(dis < sqrt(L*L))continue;

        if((objs[no].polys[nv].n*view+objs[no].polys[nv].d)*(objs[no].polys[nv].n*lights[i].loc+objs[no].polys[nv].d) > 0)
        {
            tmp = tmp + 1.0/sqrt(L*L)*L*n*objs[no].Kds;
            tmp = tmp + pow(1/sqrt((L+V)*(L+V))*(L+V)*n, 30)*objs[no].Ks;
        }
        ans = ans+tmp*lights[i].col;
    }
    return ans;
}

Color Ilocal(Point3D view, Point3D pos, int no, int nv, bool in)
{
    Color ans(0.0, 0.0, 0.0);
    for(size_t i=0;i<lights.size();i++)
    {
        int no1, nv1;
        Point3D p, L=lights[i].loc-pos, V=view-pos, tmp(0.0, 0.0, 0.0), n;
        if(L*objs[no].polys[nv].n < 0)n=-1*objs[no].polys[nv].n;else n=objs[no].polys[nv].n;
        double dis, n1=10, n2=17;
        if(!in)swap(n1, n2);
        kdtree->check(objs, pos, L, no1, nv1, p, dis);
        //CheckAll(pos, L, no1, nv1, p, dis);
        if(dis < sqrt(L*L))continue;
        V = (1.0/sqrt(V*V))*V;

        if((objs[no].polys[nv].n*view+objs[no].polys[nv].d)*(objs[no].polys[nv].n*lights[i].loc+objs[no].polys[nv].d) < 0)
        {
            if(!in && asin(n2/n1) < acos(V*n))continue;
            tmp = tmp - 1.0/sqrt(L*L)*L*n*objs[no].Kdt;
            tmp = tmp + pow((n1-n2)/fabs(n1-n2)/sqrt((n2*L+n1*V)*(n2*L+n1*V))*(n2*L+n1*V)*n, 30)*objs[no].Kt;
        }
        ans = ans+tmp*lights[i].col;
    }
    return ans;
}

Color RayTracing(Point3D view, Point3D ray, double weight, bool in)
{
    const double MinWeight = 1e-3;
    if(weight < MinWeight) return Color(0.0, 0.0, 0.0);
    double dis;
    int no=-1, nv=-1;
    Point3D p, n;
    if(!kdtree->check(objs, view, ray, no, nv, p, dis))return Color(0.0, 0.0, 0.0);
    //if(!CheckAll(view, ray, no, nv, p, dis))return Color(0.0, 0.0, 0.0);
    Point3D V = view-p;V = (1.0/sqrt(V*V))*V;
    if(V*objs[no].polys[nv].n < 0)n=-1*objs[no].polys[nv].n;else n=objs[no].polys[nv].n;

    Point3D R = 2*(V*n)*n-V;
    Point3D delta = eps/sqrt(ray*ray)*ray;
    Color ans = Ilocal(view, p-delta, no, nv)+Ilocal(view, p+delta, no, nv, in);
    Color Rcolor;
    if(objs[no].Ks1*objs[no].Ks1 > eps)Rcolor = objs[no].Ks1*RayTracing(p-delta, R, weight*0.6, in);
    ans = ans + Rcolor;
    double cos1 = n*V, np = in?(17.0/10.0):(10.0/17.0);
    if(1-1/(np*np)*(1-cos1*cos1) < eps)return ans;// = ans + 0.5*Rcolor;
    double cos2 = sqrt(1-(1/(np*np)*(1-cos1*cos1)));
    Point3D T = -1/np*V+(cos1/np-cos2)*n;
    if(objs[no].Kt1*objs[no].Kt1 > eps)ans = ans + objs[no].Kt1*RayTracing(p+delta, T, weight*0.4, !in);
    return ans;
}

int main( int argc, char** argv )
{
//    if(argc < 2)
//    {
//        printf("Usage %s inputFile\n", argv[0]);
//        return 0;
//    }
    freopen("test.txt", "r", stdin);

    Mat image;
    int n, m, sizex, sizey;
    Point3D view, sc_tmp;

    scanf("%d", &n);

    for(int i=0;i<n;i++)
    {
        objs.push_back(Object());
        string tmp;
        cin>>tmp;
        Point3D loc, rotate;
        double times;
        scanf("%lf", &times);
        scanf("%lf%lf%lf", &loc.x, &loc.y, &loc.z);
        scanf("%lf%lf%lf", &rotate.x, &rotate.y, &rotate.z);
        scanf("%lf%lf%lf", &objs[i].Kds.x, &objs[i].Kds.y, &objs[i].Kds.z);
        scanf("%lf%lf%lf", &objs[i].Ks.x, &objs[i].Ks.y, &objs[i].Ks.z);
        scanf("%lf%lf%lf", &objs[i].Ks1.x, &objs[i].Ks1.y, &objs[i].Ks1.z);
        scanf("%lf%lf%lf", &objs[i].Kdt.x, &objs[i].Kdt.y, &objs[i].Kdt.z);
        scanf("%lf%lf%lf", &objs[i].Kt.x, &objs[i].Kt.y, &objs[i].Kt.z);
        scanf("%lf%lf%lf", &objs[i].Kt1.x, &objs[i].Kt1.y, &objs[i].Kt1.z);
        if(!objs[i].readfile(tmp, times, loc, rotate)) return 0;
    }
    kdtree = new KdTree(objs, 0);

    scanf("%d", &m);
    for(int i=0;i<m;i++)
    {
        Point3D loc;
        Color col;
        scanf("%lf%lf%lf", &loc.x, &loc.y, &loc.z);
        scanf("%lf%lf%lf", &col.r, &col.g, &col.b);
        lights.push_back(Light(loc, col));
    }

    scanf("%lf%lf%lf", &view.x, &view.y, &view.z);
    vector<Point3D> sc;
    for(int i=0;i<3;i++)
    {
        scanf("%lf%lf%lf", &sc_tmp.x, &sc_tmp.y, &sc_tmp.z);
        sc.push_back(sc_tmp);
    }
    scanf("%d%d", &sizex, &sizey);

    image.create(sizex, sizey, CV_8UC3);

    Polygon screen(sc);
    Point2D Vy_2=(1.0/sizey)*screen.points2d[1], Vx_2=(1.0/sizex)*(screen.points2d[2]-screen.points2d[1]);
    Point3D Vx(Vx_2.x, Vx_2.y, 0), Vy(Vy_2.x, Vy_2.y, 0);

    for(int i=0;i<sizex;i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for(int j=0;j<sizey;j++)
        {
            Point3D ray = screen.rotate_r*(i*Vx+j*Vy) - view;
            drawPixel(image, i, j, RayTracing(view, ray, 1, false));
        }
    }
    printf("Save image to %s\n", "output.jpg");
    imwrite("output.jpg", image);
//    namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
//    imshow( "Display Image", image );
//    waitKey(0);

    return 0;
}


