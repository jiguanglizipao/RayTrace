#include <cv.h>
#if CV_VERSION_MAJOR == 3                                                                                       //支持CV2和CV3
    #include <opencv2/highgui.hpp>
#else
    #include <opencv/highgui.h>
#endif
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "point.h"
#include "polygon.h"
#include "object.h"

using namespace std;
using namespace cv;

struct Color                                                                                                    //颜色类
{
    int r, g, b;
    Color(int _r = 0, int _g = 0, int _b = 0)
        :r(_r), g(_g), b(_b)
    {
    }
};

void drawPixel(cv::Mat &image, int x, int y, Color color) {                                                     //在x,y绘制color
    image.at<cv::Vec3b>(x,y)[2] = color.r;
    image.at<cv::Vec3b>(x,y)[1] = color.g;
    image.at<cv::Vec3b>(x,y)[0] = color.b;
}

int main( int argc, char** argv )
{
    Mat image;                                                                                                  //创建读入物体
    int sizex, sizey, times;
    printf("Please input size of image & times of the obj & .obj file name.\nExample: 1080 1920 400 sphere.obj\n");
    scanf("%d%d%d", &sizex, &sizey, &times);
    string tmp;
    cin>>tmp;
    image.create(sizex, sizey, CV_8UC3);
    Object obj;
    if(!obj.readfile(tmp, times))return 0;
    #pragma omp parallel for                                                                                    //OpenMP
    for(int i=0;i<sizex;i++)
    {
        for(int j=0;j<sizey;j++)
        {
            double md = -1e20;
            drawPixel(image, i, j, Color(255, 255, 255));                                                       //每个像素先设成白色
            for(size_t k=0;k<obj.polys.size();k++)
            {
                Polygon::State state = obj.polys[k].checkInside(Point2D(i-sizex/2, j-sizey/2));                 //获得像素的状态
                if(state == Polygon::inside || state == Polygon::online)                                        //改变像素的颜色
                {
                    double depth = obj.polys[k].getDepth(Point2D(i-sizex/2, j-sizey/2));
                    if(depth < md-eps)continue;
                    md = depth;
                    if(state == Polygon::inside)
                        drawPixel(image, i, j, Color(128, 128, 128));
                    else
                        drawPixel(image, i, j, Color(0, 0, 0));
                }
                //if(state == Polygon::error)continue;
            }
        }
    }
    printf("Save image to %s", "output.jpg");
    imwrite("output.jpg", image);
    namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
    imshow( "Display Image", image );
    waitKey(0);

    return 0;
}


