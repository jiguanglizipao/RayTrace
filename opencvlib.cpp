#include <cv.hpp>
#if CV_VERSION_MAJOR == 3
#include <opencv2/highgui.hpp>
#else
#include <opencv/highgui.h>
#endif
#include "opencvlib.h"
using namespace cv;


void loadImage(unsigned char *&tex, int &sizex, int &sizey, const char file[], float times)
{
    FILE *ft = fopen(file, "r");
    if(!ft)return;else fclose(ft);
    Mat image = imread(file);
    Size size = Size(image.cols*times, image.rows*times);
    Mat im(size, CV_8UC3);
    resize(image, im, size);
    sizex = image.cols*times;
    sizey = image.rows*times;
    int nr = im.rows;
    int nc = im.cols * im.channels();
    tex = new unsigned char[nr*nc];
    unsigned char *t=tex;
    for(int i=0;i<nr;i++)
    {
        const uchar* data= im.ptr<uchar>(i);
        for(int j=0;j<nc;j++)*(t++) = data[j];
    }
}
