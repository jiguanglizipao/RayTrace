#ifndef OBJECT_H
#define OBJECT_H
#include "point.h"
#include "polygon.h"
#include <vector>
#include <string>

struct Object                                               //物体类
{
    std::vector<Polygon> polys;                             //物体的面
    std::vector<Point3D> points;                            //物体的顶点
    int times;                                              //物体的倍数
    Object(std::string fi = "", int _times = 1);
    bool readfile(std::string fi = "", int _times = 1);     //读入obj文件
};

#endif
