#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <cv.hpp>
#if CV_VERSION_MAJOR == 3
#include <opencv2/highgui.hpp>
#else
#include <opencv/highgui.h>
#endif
#include "point.h"
#include "sphere.h"
#include "polygon.h"
#include "object.h"
#include "kdtree.h"
#include "cuda_raytrace.h"

using namespace std;
using namespace cv;

vector<Object> objs;
vector<Sphere> spheres;
KdTree *kdtree;

inline double norm(double x)
{
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

void drawPixel(cv::Mat &image, int x, int y, Point3D color) {
    Point3D ans = color;
    ans.x=norm(ans.x), ans.y=norm(ans.y), ans.z=norm(ans.z);
    image.at<cv::Vec3b>(x,y)[2] = int(pow(ans.x, 1 / 2.2) * 255 + .5);
    image.at<cv::Vec3b>(x,y)[1] = int(pow(ans.y, 1 / 2.2) * 255 + .5);
    image.at<cv::Vec3b>(x,y)[0] = int(pow(ans.z, 1 / 2.2) * 255 + .5);
}

Point3D radiance(Ray r, int depth, bool into)
{
    double ts, to=1e20;
    int ids, ido, idv;
    bool fs, fo;
    fo = kdtree->check(objs, r, ido, idv, to);
    fs = sphere_intersect(&spheres[0], spheres.size(), r, ids, ts);
    if(!fo && !fs)return Point3D();
    Point3D x, n, nl, f, lig;
    RType type;
//    if(fo)printf("%lf %lf %d %d\n", ts, to, ido, idv);
    if(ts < to || !fo)
    {
        const Sphere & obj = spheres[ids];
        x = r.o + r.d * ts, n = (x - obj.pos).norm(), nl = ((n*r.d) < 0 ? n : n * -1);
        f = obj.col, type = obj.type, lig = obj.lig;
    }
    else
    {
        double u, v;
        objs[ido].polys[idv].intersect(r, u, v);
        x = r.o + r.d * to;
        n = objs[ido].polys[idv].getn(x);
        nl = ((n*r.d) < 0 ? n : n * -1);
        f = objs[ido].polys[idv].getcol(u, v);
        type = objs[ido].type, lig = objs[ido].lig;
    }

    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;	// max refl
    if (++depth > 5)
        if (drand48() < p)
            f = f * (1 / p);
        else
            return lig;


    if (type == DIFF)
    {
        double r1 = 2 * M_PI * drand48(), r2 = drand48(), r2s = sqrt(r2);
        Point3D w = nl, u = ((fabs(w.x) > .1 ? Point3D(0, 1, 0) : Point3D(1, 0, 0)) % w).norm(), v = w % u;
        Point3D d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return lig + f.mult(radiance(Ray(x-r.d*eps, d), depth, into));
    } else if (type == SPEC)	// Ideal SPECULAR reflection
        return lig + f.mult(radiance(Ray(x-r.d*eps, r.d - n * 2 * (n*r.d)), depth, into));

    Ray reflRay(x-r.d*eps, r.d - n * 2 * (n*r.d));	// Ideal dielectric REFRACTION
    double nc = 1, nt = 1.7, nnt = into ? nc / nt : nt / nc, ddn = r.d*nl, cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)	// Total internal reflection
        return lig + f.mult(radiance(reflRay, depth, into));
    Point3D tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b),   c = 1 - (into ? -ddn : (tdir*n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return lig + f.mult(depth > 2 ? (drand48() < P ? radiance(reflRay, depth, into) * RP : radiance(Ray(x+r.d*eps, tdir), depth, !into) * TP) : 
                        radiance(reflRay, depth, into) * Re + radiance(Ray(x+r.d*eps, tdir), depth, !into) * Tr);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    int myid, mpin;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &mpin);
    srand48(time(NULL)+myid);
    FILE *fi = fopen("test.txt", "r");

    int n, m, sizex, sizey, c_samps, g_samps, NUM_PER_NODE, CU_STACK_SIZE;
    fscanf(fi, "%d", &n);
    for(int i=0;i<n;i++)
    {
        objs.push_back(Object());
        char buf[256], buf2[256];
        fscanf(fi, "%s", buf);
        Point3D loc, rotate;
        double times;
        fscanf(fi, "%lf", &times);
        fscanf(fi, "%lf%lf%lf", &loc.x, &loc.y, &loc.z);
        fscanf(fi, "%lf%lf%lf", &rotate.x, &rotate.y, &rotate.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].lig.x, &objs[i].lig.y, &objs[i].lig.z);
        fscanf(fi, "%lf%lf%lf", &objs[i].col.x, &objs[i].col.y, &objs[i].col.z);
        fscanf(fi, "%s", buf2);
        if(buf2[0] == 'D')objs[i].type = DIFF;
        else if(buf2[0] == 'S')objs[i].type = SPEC;
        else objs[i].type = REFR;
        if(!objs[i].readfile(string(buf), times, loc, rotate)) return 0;
    }
    kdtree = new KdTree(objs, 0);

    fscanf(fi, "%d", &m);
    for(int i=0;i<m;i++)
    {
        double rad;
        Point3D pos, lig, col;
        RType type;
        fscanf(fi, "%lf", &rad);
        fscanf(fi, "%lf%lf%lf", &pos.x, &pos.y, &pos.z);
        fscanf(fi, "%lf%lf%lf", &lig.x, &lig.y, &lig.z);
        fscanf(fi, "%lf%lf%lf", &col.x, &col.y, &col.z);
        char buf[256];
        fscanf(fi, "%s", buf);
        if(buf[0] == 'D')type = DIFF;
        else if(buf[0] == 'S')type = SPEC;
        else type = REFR;
        spheres.push_back(Sphere(rad, pos, lig, col, type));
    }

    fscanf(fi, "%d%d%d%d%d%d", &sizex, &sizey, &c_samps, &g_samps, &NUM_PER_NODE, &CU_STACK_SIZE);
    Point3D view, dir;
    fscanf(fi, "%lf%lf%lf", &view.x, &view.y, &view.z);
    fscanf(fi, "%lf%lf%lf", &dir.x, &dir.y, &dir.z);


    //int sizex = 768, sizey = 1024, samps = 8;	// # samples
    Ray cam(view, dir.norm());	// cam pos, dir
    Point3D cy = Point3D(sizey * .5035 / sizex, 0, 0), cx = (cy % cam.d).norm() * -.5035;

    int dev_num = check_device();
    int nodeid = myid%NUM_PER_NODE;
    int dev = -1, samps=c_samps;
    //if(nodeid < dev_num)dev=nodeid, samps=g_samps;
    if(dev != -1)
    {
        print_device_info(myid, dev);
        cudaSetDevice(dev);
    #if CUDART_VERSION >= 4000
        cudaDeviceSetLimit(cudaLimitStackSize, CU_STACK_SIZE);
    #else
        cudaThreadSetLimit(cudaLimitStackSize, CU_STACK_SIZE);
    #endif
        copyToDevice(spheres, objs, kdtree);
    }

    MPI_Request req;
    Point3D* col;

    printf("myid=%d %s samps=%d\n", myid, (dev==-1)?"CPU":"GPU", samps);

    int start = (myid%2 == 0)?0:(sizex/10*7), end = (myid%2 == 0)?(sizex/10*7):sizex;
    col = new Point3D[sizey*(end-start)]();

    for(int x=start;x<end;x++)
    {
        printf("myid=%d %s calc X=%d time=%lf\n", myid, (dev==-1)?"CPU":"GPU", x, MPI_Wtime()-start_time);
        //#pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < sizey; y++)
        {	// Loop cols
            if(dev != -1)
            {
                Point3D r = cuda_radiance(x, y, sizex, sizey, samps, cx, cy, cam);
                col[(x-start)*sizey+y] = col[(x-start)*sizey+y] + r;
            }
            else
            {
                Point3D r;
                for(int sx=0;sx<2;sx++)for(int sy=0;sy<2;sy++,r=Point3D())
                {   // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * drand48(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * drand48(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Point3D d = cx * (((sx + .5 + dx) / 2 + x) / sizex - .5) + cy * (((sy + .5 + dy) / 2 + y) / sizey - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, true)*(1.0/samps);
                    }   // Camera rays are pushed ^^^^^ forward to start in interior
                    col[(x-start)*sizey+y] = col[(x-start)*sizey+y] + Point3D(norm(r.x), norm(r.y), norm(r.z)) * .25;
                }
            }
        }
    }
    
    MPI_Isend(col, 3*sizey*(end-start), MPI_DOUBLE, 0, myid, MPI_COMM_WORLD, &req);
    printf("myid=%d %s Finish\n", myid, (dev==-1)?"CPU":"GPU");

    Mat image;
    if(!myid)
    {
        image.create(sizex, sizey, CV_8UC3);
        Point3D *buff = new Point3D[sizex*sizey](), *ans = new Point3D[sizex*sizey]();
        for(int k=0;k<mpin;k++)
        {
            printf("Getdata id=%d\n", k);
            int start = (k%2 == 0)?0:(sizex/10*7), end = (k%2 == 0)?(sizex/10*7):sizex;
            MPI_Recv(buff, 3*sizey*(end-start), MPI_DOUBLE, k, k, MPI_COMM_WORLD, &status);
            for(int i=start;i<end;i++)
                for(int j=0;j<sizey;j++)
                    ans[i*sizey+j] = ans[i*sizey+j]+buff[(i-start)*sizey+j]*(2.0/mpin);
        }
        for(int i=0;i<sizex;i++)
            for(int j=0;j<sizey;j++)
                drawPixel(image, i, j, ans[i*sizey+j]);
        delete [] buff;
        delete [] ans;
        printf("Save image to %s\n", "output.bmp");
        imwrite("output.bmp", image);
        if(!myid)printf("%lfs\n", MPI_Wtime()-start_time);
    }
    MPI_Wait(&req, &status);
    MPI_Finalize();
    delete [] col;
    if(!myid)
    {
        namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
        imshow( "Display Image", image );
        waitKey(0);
    }
    return 0;

}

