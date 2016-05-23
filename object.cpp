#include "object.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
using namespace std;
using namespace cv;

bool Object::readfile(std::string filename, double _times, Point3D _loc, Point3D rotate)
{
    times = _times;
    loc = _loc;
    polys.clear();
    points.clear();
    tex.resize(1, NULL);
    cols.resize(1, col);
    name.resize(1);
    vector<Point3D> vn, vt;
    vn.clear(), vt.clear();
    vt.resize(1);
    rotate = rotate*(M_PI/180);
    Matrix t = Matrix::Rx(cos(rotate.x), sin(rotate.x))*Matrix::Ry(cos(rotate.y), sin(rotate.y))*Matrix::Rz(cos(rotate.z), sin(rotate.z));
    FILE *fp = fopen(filename.c_str(), "r");
    if(!fp){
        printf("Can't open %s\n", filename.c_str());
        return false;
    }
    char buf[512];
    char str[512];
    int lines = 0;
    vector<vector<int> >face, N, T;
    vector<int>Tex;
    int ntex=0;
    while(fscanf(fp, "%s", buf) != EOF)
    {
        lines++;
        switch(buf[0])
        {
        case '#':				/* comment */
            /* eat up rest of line */
            fgets(buf, sizeof(buf), fp);
            break;
        case 'u':
            fgets(buf, sizeof(buf), fp);
            sscanf(buf, "%s", str);
            ntex = 0;
            for(int i=0;i<tex.size();i++)if(name[i] == string(str))ntex = i;
            break;

        case 'm':				/* comment */
            /* eat up rest of line */
            fgets(buf, sizeof(buf), fp);
            sscanf(buf, "%s", str);
            FILE *fm = fopen(str, "r");
            while(fscanf(fm, "%s", buf) != EOF)
            {
                switch(buf[0])
                {
                case 'n':
                    fgets(buf, sizeof(buf), fm);
                    sscanf(buf, "%s", str);
                    name.push_back(string(str));
                    tex.push_back(NULL);
                    cols.push_back(col);
                    break;
                case 'K':
                    switch(buf[1])
                    {
                    case 'a':
                        fgets(buf, sizeof(buf), fm);
                        break;
                    case 'd':
                        if(type == DIFF)
                        {
                            Point3D v;
                            if(fscanf(fm, "%lf%lf%lf",&v.x, &v.y, &v.z)==3)
                            {
                            }
                            else
                            {
                                printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lines);
                                return false;
                            }
                            cols.back()=v;
                        }
                        else fgets(buf, sizeof(buf), fm);
                        break;
                    case 's':
                        if(type != DIFF)
                        {
                            Point3D v;
                            if(fscanf(fm, "%lf%lf%lf",&v.x, &v.y, &v.z)==3)
                            {
                            }
                            else
                            {
                                printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lines);
                                return false;
                            }
                            cols.back()=v;
                        }
                        else fgets(buf, sizeof(buf), fm);
                        break;
                    default:
                        fgets(buf, sizeof(buf), fm);
                        break;
                    }
                    break;
                case 'm':
                    if((!strcmp(buf, "map_Kd") && type == DIFF) || (!strcmp(buf, "map_Ks") && type != DIFF))
                    {
                        memset(buf, 0, sizeof(buf));
                        fgets(buf, sizeof(buf), fm);
                        sscanf(buf, "%s", str);
                        FILE *ft = fopen(str, "r");
                        if(!ft)break;else fclose(ft);
                        Mat image = imread(str);
                        Size size = Size(image.cols*times,image.rows*times);
                        tex.back() = new Mat(size, CV_8UC3);
                        resize(image, *tex.back(), size);
                    }
                    else fgets(buf, sizeof(buf), fm);
                    break;
                default:
                    fgets(buf, sizeof(buf), fm);
                    break;
                }
            }
            fclose(fm);
            break;
        case 'v':				/* v, vn, vt */
            switch(buf[1])
            {
            case '\0':			    /* vertex */
            {
                Point3D v;
                if(fscanf(fp, "%lf%lf%lf",&v.x, &v.y, &v.z)==3)
                {
                    v.x*=times, v.y*=times, v.z*=times;
                    v = t*v;
                    points.push_back(v+loc);
                }
                else
                {
                    printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lines);
                    return false;
                }
                break;
            }
            case 'n':
            {
                Point3D v;
                if(fscanf(fp, "%lf%lf%lf",&v.x, &v.y, &v.z)==3)
                {
                    v = t*v;
                    vn.push_back(v.norm());
                }
                else
                {
                    printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lines);
                    return false;
                }
                break;
            }
            case 't':
            {
                Point3D v;
                if(fscanf(fp, "%lf%lf",&v.x, &v.y)==2)
                {
                    vt.push_back(v);
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

        case 'f':				/* face */
        {
            vector<int>t, tmp, n;
            t.resize(3);tmp.resize(3);n.resize(3);
            if(fscanf(fp, "%s", str)!=1)
            {
                printf("Error: Wrong Face at Line %d\n",lines);
                return false;
            }
            if (strstr(str, "//"))
            {
                /* v//n */
                if( sscanf(str, "%d//%d", &tmp[0], &n[0]) ==2  &&
                    fscanf(fp, "%d//%d", &tmp[1], &n[1]) ==2  &&
                    fscanf(fp, "%d//%d", &tmp[2], &n[2]) ==2)
                {
                }
                else
                {
                    printf("Error: Wrong Face at Line %d\n",lines);
                    return false;
                }

            }
            else if (sscanf(str, "%d/%d/%d", &tmp[0], &t[0], &n[0]) == 3)
            {
                /* v/t/n */
                if( fscanf(fp, "%d/%d/%d", &tmp[1], &t[1], &n[1]) ==3 &&
                    fscanf(fp, "%d/%d/%d", &tmp[2], &t[2], &n[2]) ==3 )
                {
                }
                else
                {
                    printf("Error: Wrong Face at Line %d\n",lines);
                    return false;
                }
            }
            else if (sscanf(str, "%d/%d", &tmp[0], &t[0]) == 2)
            {
                /* v/t */
                if( fscanf(fp, "%d/%d", &tmp[1], &t[1]) ==2 &&
                    fscanf(fp, "%d/%d", &tmp[2], &t[2]) ==2 )
                {
                }
                else
                {
                    printf("Error: Wrong Face at Line %d\n",lines);
                    return false;
                }
            }
            else
            {
                /* v */
                if( sscanf(str, "%d", &tmp[0]) ==1 &&
                    fscanf(fp, "%d", &tmp[1])  ==1 &&
                    fscanf(fp, "%d", &tmp[2])  ==1 )
                {
                }
                else
                {
                    printf("Error: Wrong Face at Line %d\n",lines);
                    return false;
                }
            }
            tmp[0]--;tmp[1]--;tmp[2]--;
            n[0]--;n[1]--;n[2]--;
            face.push_back(tmp);
            N.push_back(n);
            T.push_back(t);
            Tex.push_back(ntex);
            break;
        }
        default:
            /* eat up rest of line */
            fgets(buf, sizeof(buf), fp);
            break;
        }
    }

    for(std::size_t i=0;i<face.size();i++)
    {
        Point3D tmp[3], n[3], t[3];
        int ntex;
        for(std::size_t j=0;j<3;j++)
        {
            if(face[i][j] >= points.size() || face[i][j] < 0)
            {
                printf("Error face.\n");
                points.clear();
                polys.clear();
                return false;
            }
            tmp[j] = points[face[i][j]];
            n[j] = vn[N[i][j]];
            t[j] = vt[T[i][j]];
            ntex = Tex[i];
        }
        polys.push_back(Polygon(tmp, n, t, lig, cols[ntex], tex[ntex]));
    }
    return true;
}

