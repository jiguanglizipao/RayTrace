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
//#include "polygon.h"
//#include "obj.colt.h"
//#include "kdtree.h"
//#include "light.h"

using namespace std;
using namespace cv;

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

enum RType { DIFF, SPEC, REFR };	// material types, used in radiance()
struct Sphere
{
    double rad;
    Point3D pos, lig, col;
    RType type;
    Sphere(double rad_, Point3D pos_, Point3D lig_, Point3D col_, RType type_)
        :rad(rad_), pos(pos_), lig(lig_), col(col_), type(type_) {
    }
    double intersect(const Ray & r) const
    {
        Point3D op = pos - r.o;
        double t, b = op * r.d, det = b * b - op * op + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};
Sphere spheres[] = {		//Scene: radius, position, emission, color, material
                            Sphere(1e5, Point3D(1e5 + 1, 40.8, 81.6), Point3D(), Point3D(.75, .25, .25), DIFF),	//Left
                            Sphere(1e5, Point3D(-1e5 + 99, 40.8, 81.6), Point3D(), Point3D(.25, .25, .75), DIFF),	//Rght
                            Sphere(1e5, Point3D(50, 40.8, 1e5), Point3D(), Point3D(.75, .75, .75), DIFF),	//Back
                            Sphere(1e5, Point3D(50, 40.8, -1e5 + 170), Point3D(), Point3D(.25, .75, .25), DIFF),	//Frnt
                            Sphere(1e5, Point3D(50, 1e5, 81.6), Point3D(), Point3D(.75, .75, .75), DIFF),	//Botm
                            Sphere(1e5, Point3D(50, -1e5 + 81.6, 81.6), Point3D(), Point3D(.75, .75, .75), DIFF),	//Top
                            Sphere(16.5, Point3D(27, 16.5, 47), Point3D(), Point3D(1, 1, 1) * .999, SPEC),	//Mirr
                            Sphere(16.5, Point3D(73, 16.5, 78), Point3D(), Point3D(1, 1, 1) * .999, REFR),	//Glas
                            Sphere(600, Point3D(50, 681.6 - .27, 81.6), Point3D(3, 3, 3), Point3D(), DIFF)	//Lite
                   };

inline bool intersect(const Ray & r, double &t, int &id)
{
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int (n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}

Point3D radiance(const Ray & r, int depth, bool into)
{
    double t;		// distance to intersection
    int id = 0;		// id of intersected obj.colt
    if (!intersect(r, t, id))return Point3D();	// if miss, return black

    const Sphere & obj = spheres[id];	// the hit obj.colt
    Point3D x = r.o + r.d * t, n = (x - obj.pos).norm(), nl = (n*r.d) < 0 ? n : n * -1, f = obj.col;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;	// max refl
    if (++depth > 5)
        if (drand48() < p)
            f = f * (1 / p);
        else
            return obj.lig;

    if (obj.type == DIFF)
    {
        double r1 = 2 * M_PI * drand48(), r2 = drand48(), r2s = sqrt(r2);
        Point3D w = nl, u = ((fabs(w.x) > .1 ? Point3D(0, 1, 0) : Point3D(1, 0, 0)) % w).norm(), v = w % u;
        Point3D d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.lig + f.mult(radiance(Ray(x, d), depth, into));
    } else if (obj.type == SPEC)	// Ideal SPECULAR reflection
        return obj.lig + f.mult(radiance(Ray(x, r.d - n * 2 * (n*r.d)), depth, into));

    Ray reflRay(x, r.d - n * 2 * (n*r.d));	// Ideal dielectric REFRACTION
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d*nl, cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)	// Total internal reflection
        return obj.lig + f.mult(radiance(reflRay, depth, into));
    Point3D tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : (tdir*n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj.lig + f.mult(depth > 2 ? (drand48() < P ? radiance(reflRay, depth, into) * RP : radiance(Ray(x, tdir), depth, !into) * TP) :
                                        radiance(reflRay, depth, into) * Re + radiance(Ray(x, tdir), depth, !into) * Tr);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    double start = MPI_Wtime();
    int myid, mpin;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &mpin);

    int sizex = 768, sizey = 1024, samps = 8;	// # samples
    Ray cam(Point3D(50, 52, 295.6), Point3D(0, -0.042612, -1).norm());	// cam pos, dir
    Point3D cy = Point3D(sizey * .5035 / sizex, 0, 0), cx = (cy % cam.d).norm() * -.5035;

    int *spl = new int[mpin+1], per=sizex/mpin, res=sizex-mpin*per;
    memset(spl, 0, sizeof(spl));
    for(int i=mpin-1;i>=0;i--,res--)spl[i+1]=per+(res>0);
    for(int i=1;i<=mpin;i++)spl[i]=spl[i-1]+spl[i];

    vector<MPI_Request*> req;
    vector<Point3D*> col;

    printf("myid=%d start=%d end=%d\n", myid, spl[myid], spl[myid+1]);

    for(int x=spl[myid];x<spl[myid+1];x++)
    {
        printf("%d\n", x);
        col.push_back(new Point3D[sizey]);
        req.push_back(new MPI_Request[3*sizey]);
        #pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < sizey; y++)
        {	// Loop cols
            Point3D r = Point3D(0);
            for (int sy = 0; sy < 2; sy++)	// 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++)
                {	// 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * drand48(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * drand48(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Point3D d = cx * (((sx + .5 + dx) / 2 + x) / sizex - .5) + cy * (((sy + .5 + dy) / 2 + y) / sizey - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, true) * (1. / samps);
                    }	// Camera rays are pushed ^^^^^ forward to start in interior
                    col.back()[y] = col.back()[y] + Point3D(norm(r.x), norm(r.y), norm(r.z)) * .25;
                }
        }
    }

    for(int i=spl[myid];i<spl[myid+1];i++)
    {
        for(int j=0;j<sizey;j++)
        {
            MPI_Isend(&col[i-spl[myid]][j].x, 1, MPI_DOUBLE, 0, (i*sizey+j)*3+0, MPI_COMM_WORLD, &req[i-spl[myid]][3*j+0]);
            MPI_Isend(&col[i-spl[myid]][j].y, 1, MPI_DOUBLE, 0, (i*sizey+j)*3+1, MPI_COMM_WORLD, &req[i-spl[myid]][3*j+1]);
            MPI_Isend(&col[i-spl[myid]][j].z, 1, MPI_DOUBLE, 0, (i*sizey+j)*3+2, MPI_COMM_WORLD, &req[i-spl[myid]][3*j+2]);
        }
    }
    printf("myid=%d Finish\n", myid);

    Mat image;
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
                    Point3D tmp(0,0,0);
                    MPI_Recv(&(tmp.x), 1, MPI_DOUBLE, k, (i*sizey+j)*3+0, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(tmp.y), 1, MPI_DOUBLE, k, (i*sizey+j)*3+1, MPI_COMM_WORLD, &status);
                    MPI_Recv(&(tmp.z), 1, MPI_DOUBLE, k, (i*sizey+j)*3+2, MPI_COMM_WORLD, &status);
                    drawPixel(image, i, j, tmp);
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

