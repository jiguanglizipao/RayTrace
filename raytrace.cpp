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
#include "light.h"

using namespace std;
using namespace cv;

vector<Object> objs;
vector<Light> lights;

void drawPixel(cv::Mat &image, int x, int y, Color color) {
    Color ans = color;
    if(ans.b > 1 || ans.r > 1 || ans.g > 1)
        printf("%f %f %f\n", ans.r, ans.g, ans.b);
    image.at<cv::Vec3b>(x,y)[2] = int(color.r*255);
    image.at<cv::Vec3b>(x,y)[1] = int(color.g*255);
    image.at<cv::Vec3b>(x,y)[0] = int(color.b*255);
}

bool CheckAll(Point3D view, Point3D ray, int &no, int &nv, Point3D &p, float &dis)
{
    dis=1.1e10;
    for(size_t i=0;i<objs.size();i++)
    {
        for(size_t j=0;j<objs[i].polys.size();j++)
        {
            float k = -(objs[i].polys[j].d+objs[i].polys[j].n*view)/(objs[i].polys[j].n*ray) - eps;
            if(k < eps)continue;
            Point3D tmp = view+k*ray;
            if(objs[i].polys[j].checkInside(tmp) == Polygon::inside)
                if(k < dis)dis=k, no=i, nv=j, p=tmp;
        }
    }
    return dis < 1e10;
}

Color Ilocal(Point3D view, Point3D pos, int no, int nv, bool in)
{
    Color ans(0.0, 0.0, 0.0);
    for(size_t i=0;i<lights.size();i++)
    {
        int no1, nv1;
        Point3D p, L=lights[i].loc-pos, V=view-pos, tmp(0.0, 0.0, 0.0);
        float dis, n1=3, n2=5;
        if(in)swap(n1, n2);
        if(CheckAll(pos, L, no1, nv1, p, dis))continue;

        if((objs[no].polys[nv].n*view+objs[no].polys[nv].d)*(objs[no].polys[nv].n*lights[i].loc+objs[no].polys[nv].d) > 0)
        {
            tmp = tmp + 0.5/sqrt(L*L)*L*objs[no].polys[nv].n*objs[no].Ks;
            tmp = tmp + 0.5*pow(1/sqrt((L+V)*(L+V))*(L+V)*objs[no].polys[nv].n, 3)*objs[no].Ks;
        }
        else
        {
            tmp = tmp - 0.5/sqrt(L*L)*L*objs[no].polys[nv].n*objs[no].Kt;
            tmp = tmp + 0.5*pow((n1-n2)/fabs(n1-n2)/sqrt((n2*L+n1*V)*(n2*L+n1*V))*(n2*L+n1*V)*objs[no].polys[nv].n, 3)*objs[no].Kt;
        }
        ans = ans+tmp*lights[i].col;
    }
    return ans;
}

Color RayTracing(Point3D view, Point3D ray, float weight, bool in)
{
    const float MinWeight = 1e-3;
    if(weight < MinWeight) return Color(0.0, 0.0, 0.0);
    float dis;
    int no=-1, nv=-1;
    Point3D p;
    if(!CheckAll(view, ray, no, nv, p, dis))return Color(0.0, 0.0, 0.0);
    Point3D V = view-p;
    float thetai = acos(V*objs[no].polys[nv].n/sqrt(V*V))*2;
    Point3D R = Matrix::RotateMatrix(p, V^objs[no].polys[nv].n, thetai)*view-p;
    Point3D T = Matrix::RotateMatrix(p, V^objs[no].polys[nv].n, thetai+M_PI-asin(sin(thetai)*(in?(5/3):(3/5))))*view-p;

    return Ilocal(view, p, no, nv, in)+objs[no].Ks*RayTracing(p, R, weight*0.4, in)+objs[no].Kt*RayTracing(p, T, weight*0.3, !in);
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
        Point3D loc;
        float times;
        scanf("%f", &times);
        scanf("%f%f%f", &loc.x, &loc.y, &loc.z);
        scanf("%f%f%f", &objs[i].Ks.x, &objs[i].Ks.y, &objs[i].Ks.z);
        scanf("%f%f%f", &objs[i].Kt.x, &objs[i].Kt.y, &objs[i].Kt.z);
        if(!objs[i].readfile(tmp, times, loc)) return 0;
    }

    scanf("%d", &m);
    for(int i=0;i<m;i++)
    {
        Point3D loc;
        Color col;
        scanf("%f%f%f", &loc.x, &loc.y, &loc.z);
        scanf("%f%f%f", &col.r, &col.g, &col.b);
        lights.push_back(Light(loc, col));
    }

    scanf("%f%f%f", &view.x, &view.y, &view.z);
    vector<Point3D> sc;
    for(int i=0;i<3;i++)
    {
        scanf("%f%f%f", &sc_tmp.x, &sc_tmp.y, &sc_tmp.z);
        sc.push_back(sc_tmp);
    }
    scanf("%d%d", &sizex, &sizey);

    image.create(sizex, sizey, CV_8UC3);

    Polygon screen(sc);
    Point2D Vy_2=(1.0/sizey)*screen.points2d[1], Vx_2=(1.0/sizex)*(screen.points2d[2]-screen.points2d[1]);
    Point3D Vx(Vx_2.x, Vx_2.y, 0), Vy(Vy_2.x, Vy_2.y, 0);

    for(int i=0;i<sizex ;i++)
    {
        #pragma omp parallel for schedule(dynamic)
        for(int j=0;j<sizey;j++)
        {
            Point3D ray = screen.rotate_r*(i*Vx+j*Vy) - view;
            drawPixel(image, i, j, RayTracing(view, ray, 1, false));
        }
    }
    printf("Save image to %s", "output.jpg");
    imwrite("output.jpg", image);
    namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
    imshow( "Display Image", image );
    waitKey(0);

    return 0;
}


