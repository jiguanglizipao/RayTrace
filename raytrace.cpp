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
#include <mpi.h>
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
int lnum;

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

Color Ilocal(Point3D view, Point3D pos, int no, int nv, bool f)
{
    int lnumt = f?lnum:lnum;
    Color ans(0.0, 0.0, 0.0);
    for(size_t i=0;i<lights.size();i++)
    {
        int no1, nv1, num=0;
        Point3D p, L, V=view-pos, tmp(0.0, 0.0, 0.0), n, sum(0,0,0);
        for(int j=0;j<lnumt;j++)
        {
            double x=1e-6*(rand()%1000000), y=1e-6*(rand()%1000000), z=1e-6*(rand()%1000000), dis;
            x=lights[i].mi[0]+x*(lights[i].ma[0]-lights[i].mi[0]);
            y=lights[i].mi[1]+y*(lights[i].ma[1]-lights[i].mi[1]);
            z=lights[i].mi[2]+z*(lights[i].ma[2]-lights[i].mi[2]);
            L=Point3D(x,y,z)-pos;
            kdtree->check(objs, pos, L, no1, nv1, p, dis);
            if(dis < sqrt(L*L))continue;
            sum = sum + Point3D(x,y,z);num++;
        }
        //CheckAll(pos, L, no1, nv1, p, dis);
        if(!num)continue;
        L = 1.0/num*sum;
        if((L-pos)*objs[no].polys[nv].n < 0)n=-1*objs[no].polys[nv].n;else n=objs[no].polys[nv].n;
        if((objs[no].polys[nv].n*view+objs[no].polys[nv].d)*(objs[no].polys[nv].n*L+objs[no].polys[nv].d) > 0)
        {
            L = L-pos;
            tmp = tmp + 1.0/sqrt(L*L)*L*n*objs[no].Kds;
            tmp = tmp + pow(1/sqrt((L+V)*(L+V))*(L+V)*n, 30)*objs[no].Ks;
        }
        ans = ans+1.0*num/lnumt*tmp*lights[i].col;
    }
    return ans;
}

Color Ilocal(Point3D view, Point3D pos, int no, int nv, bool in, bool f)
{
    int lnumt = f?lnum:lnum;
    Color ans(0.0, 0.0, 0.0);
    for(size_t i=0;i<lights.size();i++)
    {
        int no1, nv1, num=0;
        Point3D p, L, V=view-pos, tmp(0.0, 0.0, 0.0), n, sum(0,0,0);
        if(L*objs[no].polys[nv].n < 0)n=-1*objs[no].polys[nv].n;else n=objs[no].polys[nv].n;
        double dis, n1=10, n2=17;
        if(!in)swap(n1, n2);
        for(int j=0;j<lnumt;j++)
        {
            double x=1e-6*(rand()%1000000), y=1e-6*(rand()%1000000), z=1e-6*(rand()%1000000), dis;
            x=lights[i].mi[0]+x*(lights[i].ma[0]-lights[i].mi[0]);
            y=lights[i].mi[1]+y*(lights[i].ma[1]-lights[i].mi[1]);
            z=lights[i].mi[2]+z*(lights[i].ma[2]-lights[i].mi[2]);
            L=Point3D(x,y,z)-pos;
            kdtree->check(objs, pos, L, no1, nv1, p, dis);
            if(dis < sqrt(L*L))continue;
            sum = sum + Point3D(x,y,z);num++;
        }
        //CheckAll(pos, L, no1, nv1, p, dis);
        L = 1.0/num*sum;
        V = (1.0/sqrt(V*V))*V;

        if((objs[no].polys[nv].n*view+objs[no].polys[nv].d)*(objs[no].polys[nv].n*L+objs[no].polys[nv].d) < 0)
        {
            L = L-pos;
            if(!in && asin(n2/n1) < acos(V*n))continue;
            tmp = tmp - 1.0/sqrt(L*L)*L*n*objs[no].Kdt;
            tmp = tmp + pow((n1-n2)/fabs(n1-n2)/sqrt((n2*L+n1*V)*(n2*L+n1*V))*(n2*L+n1*V)*n, 30)*objs[no].Kt;
        }
        ans = ans+1.0*num/lnumt*tmp*lights[i].col;
    }
    return ans;
}


Color RayTracing(Point3D view, Point3D ray, Point3D weight, bool in, int snum, bool f)
{
    const double MinWeight = 1e-3;
    double dis=1e10;
    int no=-1, nv=-1;
    Point3D p, n;
    if(!kdtree->check(objs, view, ray, no, nv, p, dis))return Color(0.0, 0.0, 0.0);
    //if(!CheckAll(view, ray, no, nv, p, dis))return Color(0.0, 0.0, 0.0);
    for(int i=0;i<lights.size();i++)if(lights[i].check_aabb(view, ray) < 1e10)return weight*lights[i].col;
    if(weight*weight < MinWeight*MinWeight) return Color(0.0, 0.0, 0.0);
    Point3D V = view-p;V = (1.0/sqrt(V*V))*V;
    if(V*objs[no].polys[nv].n < 0)n=-1*objs[no].polys[nv].n;else n=objs[no].polys[nv].n;

    Point3D R = 2*(V*n)*n-V;
    Point3D delta = eps/sqrt(ray*ray)*ray;
    Color ans = Ilocal(view, p-delta, no, nv, !f)+Ilocal(view, p+delta, no, nv, in, !f);

    if(objs[no].Ks1*objs[no].Ks1 > eps)ans = ans + RayTracing(p-delta, R, 0.6*weight*objs[no].Ks1, in, snum, f);
    if(objs[no].Ks2*objs[no].Ks2 > eps && !f)
    {
        for(int i=0;i<snum;i++)
        {
            Point3D K(1e-6*(rand()%2000000)-1, 1e-6*(rand()%2000000)-1, 1e-6*(rand()%2000000)-1);
            if((K*objs[no].polys[nv].n)*(R*objs[no].polys[nv].n) < -eps)K = -1.0*K;
            ans = ans + .1/MinWeight/snum*objs[no].Ks2*RayTracing(p-delta, K, 10*MinWeight*weight, in, snum, true);
        }
    }
    double cos1 = n*V, np = in?(17.0/10.0):(10.0/17.0);
    if(1-1/(np*np)*(1-cos1*cos1) < eps)return weight*ans;// = ans + 0.5*Rcolor;
    double cos2 = sqrt(1-(1/(np*np)*(1-cos1*cos1)));
    Point3D T = -1/np*V+(cos1/np-cos2)*n;
    if(objs[no].Kt1*objs[no].Kt1 > eps)ans = ans + RayTracing(p+delta, T, 0.4*weight*objs[no].Kt1, !in, snum, f);
    if(objs[no].Kt2*objs[no].Kt2 > eps && !f)
    {
        for(int i=0;i<snum;i++)
        {
            Point3D K(1e-6*(rand()%2000000)-1, 1e-6*(rand()%2000000)-1, 1e-6*(rand()%2000000)-1);
            if((K*objs[no].polys[nv].n)*(T*objs[no].polys[nv].n) < -eps)K = -1.0*K;
            ans = ans + .2/MinWeight/snum*objs[no].Kt2*RayTracing(p+delta, K, 10*MinWeight*weight, !in, snum, true);
        }
    }
    return weight*ans;
}

int main( int argc, char** argv )
{
    MPI_Init(&argc, &argv);
    double start = MPI_Wtime();
    int myid, mpin;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &mpin);
    //if(argc < 2)
    //{
    //    if(!myid)printf("Usage %s inputFile\n", argv[0]);
    //    return 0;
    //}

    FILE *fi = fopen("test.txt", "r");

    Mat image;
    int n, m, sizex, sizey, vnum, snum;
    Point3D view, sc_tmp;

    fscanf(fi, "%d", &n);
    for(int i=0;i<n;i++)
    {
        objs.push_back(Object());
        char buf[256];
        fscanf(fi, "%s", buf);
        Point3D loc, rotate;
        double times;
        fscanf(fi, "%lf", &times);
        fscanf(fi, "%lf%lf%lf", &loc.x, &loc.y, &loc.z);
        fscanf(fi, "%lf%lf%lf", &rotate.x, &rotate.y, &rotate.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Kds.x, &objs[i].Kds.y, &objs[i].Kds.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Ks.x, &objs[i].Ks.y, &objs[i].Ks.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Ks1.x, &objs[i].Ks1.y, &objs[i].Ks1.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Ks2.x, &objs[i].Ks2.y, &objs[i].Ks2.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Kdt.x, &objs[i].Kdt.y, &objs[i].Kdt.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Kt.x, &objs[i].Kt.y, &objs[i].Kt.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Kt1.x, &objs[i].Kt1.y, &objs[i].Kt1.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].Kt2.x, &objs[i].Kt2.y, &objs[i].Kt2.z);
        if(!objs[i].readfile(string(buf), times, loc, rotate)) return 0;
    }
    kdtree = new KdTree(objs, 0);

    fscanf(fi, "%d", &m);
    for(int i=0;i<m;i++)
    {
        Color col;
        fscanf(fi, "%lf%lf%lf", &col.r, &col.g, &col.b);
        lights.push_back(Light(col));
        fscanf(fi, "%lf%lf%lf", &lights.back().mi[0], &lights.back().mi[1], &lights.back().mi[2]);
        fscanf(fi, "%lf%lf%lf", &lights.back().ma[0], &lights.back().ma[1], &lights.back().ma[2]);
    }

    fscanf(fi, "%lf%lf%lf", &view.x, &view.y, &view.z);
    vector<Point3D> sc;
    for(int i=0;i<3;i++)
    {
        fscanf(fi, "%lf%lf%lf", &sc_tmp.x, &sc_tmp.y, &sc_tmp.z);
        sc.push_back(sc_tmp);
    }
    fscanf(fi, "%d%d%d%d%d", &sizex, &sizey, &vnum, &snum, &lnum);vnum=1;
    fclose(fi);
    Polygon screen(sc);
    Point2D Vy_2=(1.0/sizey)*screen.points2d[1], Vx_2=(1.0/sizex)*(screen.points2d[2]-screen.points2d[1]);
    Point3D Vx(Vx_2.x, Vx_2.y, 0), Vy(Vy_2.x, Vy_2.y, 0);

    int *spl = new int[mpin+1], per=sizex/mpin, res=sizex-mpin*per;
    memset(spl, 0, sizeof(spl));
    for(int i=mpin-1;i>=0;i--,res--)spl[i+1]=per+(res>0);
    for(int i=1;i<=mpin;i++)spl[i]=spl[i-1]+spl[i];

    vector<MPI_Request*> req;
    vector<Color*> col;

    printf("myid=%d start=%d end=%d\n", myid, spl[myid], spl[myid+1]);

    for(int i=spl[myid];i<spl[myid+1];i++)
    {
        printf("%d\n", i);
        col.push_back(new Color[sizey]);
        req.push_back(new MPI_Request[3*sizey]);
        #pragma omp parallel for schedule(dynamic)
        for(int j=0;j<sizey;j++)
        {
            col.back()[j]=Color(0,0,0);
            for(int k=0;k<vnum;k++)
            {
                double dx=double(rand()%100000)/100000-0.5, dy=double(rand()%100000)/100000-0.5;
                dx=0, dy=0;
                Point3D ray = screen.rotate_r*((i+dx)*Vx+(j+dy)*Vy) - view;
                col.back()[j] = col.back()[j] + RayTracing(view, ray, Point3D(1,1,1), false, snum, false);
            }
        }
    }
    for(int i=spl[myid];i<spl[myid+1];i++)
    {
        for(int j=0;j<sizey;j++)
        {
            MPI_Isend(&col[i-spl[myid]][j].r, 1, MPI_DOUBLE, 0, (i*sizey+j)*3+0, MPI_COMM_WORLD, &req[i-spl[myid]][3*j+0]);
            MPI_Isend(&col[i-spl[myid]][j].g, 1, MPI_DOUBLE, 0, (i*sizey+j)*3+1, MPI_COMM_WORLD, &req[i-spl[myid]][3*j+1]);
            MPI_Isend(&col[i-spl[myid]][j].b, 1, MPI_DOUBLE, 0, (i*sizey+j)*3+2, MPI_COMM_WORLD, &req[i-spl[myid]][3*j+2]);
        }
    }
    printf("myid=%d Finish\n", myid);
    if(!myid)
    {
        image.create(sizex, sizey, CV_8UC3);
        for(int k=0;k<mpin;k++)
        {
            printf("Getdata id=%d\n", k);
            for(int i=spl[k];i<spl[k+1];i++)
            {
                for(int j=0;j<sizey;j++)
                {
                    Color tmp(0,0,0);
                    MPI_Recv(&(tmp.r), 1, MPI_DOUBLE, k, (i*sizey+j)*3+0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(tmp.g), 1, MPI_DOUBLE, k, (i*sizey+j)*3+1, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(tmp.b), 1, MPI_DOUBLE, k, (i*sizey+j)*3+2, MPI_COMM_WORLD, &status);
                    drawPixel(image, i, j, 1.0/vnum*tmp);
                }
            }
        }
        printf("Save image to %s\n", "output.jpg");
        imwrite("output.jpg", image);
        if(!myid)printf("%lfs\n", MPI_Wtime()-start);
    }
    for(int i=0;i<req.size();i++)
        for(int j=0;j<3*sizey;j++)
            MPI_Wait(&req[i][j], &status);
     MPI_Finalize(); 
    for(int i=0;i<req.size();i++)
    {
       delete [] req[i];
       delete [] col[i];
    }
    delete [] spl;
    if(!myid)
    {
        namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
        imshow( "Display Image", image );
        waitKey(0);
    }
    return 0;
}


