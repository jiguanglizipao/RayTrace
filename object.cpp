#include "object.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>

Object::Object(std::string fi, int _times)
    :times(_times)
{
    if(fi == "")return;
    readfile(fi);
}

bool Object::readfile(std::string filename, int _times)
{
    times = _times;
    polys.clear();
    points.clear();
    FILE *fp = fopen(filename.c_str(), "r");                                                        //打开文件
    if(!fp){                                                                                        //判断是否存在
        printf("Can't open %s", filename.c_str());
        return false;
    }
    char buf[512];
    char str[512];
    int lines = 0;
    std::vector<std::vector<int> >face;                                                             //物体的面
    while(fscanf(fp, "%s", buf) != EOF)
    {
        lines++;
        switch(buf[0])
        {
        case '#':				/* comment */                                                       //注释
            /* eat up rest of line */
            fgets(buf, sizeof(buf), fp);
            break;
        case 'v':				/* v, vn, vt */                                                     //顶点
            switch(buf[1])
            {
            case '\0':			    /* vertex */
            {
                Point3D v;
                if(fscanf(fp, "%lf %lf %lf",&v.x, &v.y, &v.z)==3)
                {
                    v.x*=times, v.y*=times, v.z*=times;
                    points.push_back(v);
                }
                else
                {
                    printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lines);
                    return false;
                }
                break;
            }
            default:
                /* eat up rest of line */
                fgets(buf, sizeof(buf), fp);
                break;
            }
            break;

        case 'f':				/* face */                                                          //面
        {
            std::vector<int> tmp;
            tmp.clear();
            fgets(str, 512, fp);
            int len = 0, x;
            while(sscanf(str+len, "%d", &x)!=EOF)
            {
                tmp.push_back(x-1);
                sprintf(buf, "%d", x);
                len+=strlen(buf)+1;
            }
            face.push_back(tmp);
            break;
        }
        default:
            /* eat up rest of line */
            fgets(buf, sizeof(buf), fp);
            break;
        }
    }

    for(std::size_t i=0;i<face.size();i++)                                                          //将面与点信息写入类内
    {
        if(face[i].size() < 3)
        {
            printf("Error face.\n");
            points.clear();
            polys.clear();
            return false;
        }
        std::vector<Point3D>tmp;
        tmp.clear();
        for(std::size_t j=0;j<face[i].size();j++)
        {
            if(face[i][j] >= points.size() || face[i][j] < 0)
            {
                printf("Error face.\n");
                points.clear();
                polys.clear();
                return false;
            }
            tmp.push_back(points[face[i][j]]);
        }
        polys.push_back(Polygon(tmp));
    }
    return true;
}

