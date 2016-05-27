#ifndef KDTREE_H
#define KDTREE_H

#include <cuda_runtime.h>
#include "point.h"
#include "polygon.h"
#include "object.h"
#include <algorithm>

struct KdTreeNode
{
    int no, nv;
    float mi[3], ma[3];
    __device__ __host__ KdTreeNode()
    {
        no=nv=-1;
    }

    __device__ __host__ float min(float a, float b)
    {
        if(a < b)return a;else return b;
    }
    
    __device__ __host__ float max(float a, float b)
    {
        if(a > b)return a;else return b;
    }

    __device__ __host__ KdTreeNode(const Polygon &a, int _no, int _nv)
        :no(_no), nv(_nv)
    {
        init();
        for(size_t i=0;i<3;i++)update(a.points3d[i]);
    }

    __device__ __host__ void init(float min=-1e10, float max=1e10)
    {
        mi[0]=mi[1]=mi[2]=max;
        ma[0]=ma[1]=ma[2]=min;
    }

    __device__ __host__ void update(const Point3D &a)
    {
        mi[0] = min(mi[0], a.x);
        mi[1] = min(mi[1], a.y);
        mi[2] = min(mi[2], a.z);
        ma[0] = max(ma[0], a.x);
        ma[1] = max(ma[1], a.y);
        ma[2] = max(ma[2], a.z);
    }

    __device__ __host__ void update(const KdTreeNode &a)
    {
        for(int i=0;i<3;i++)mi[i] = min(mi[i], a.mi[i]), ma[i] = max(ma[i], a.ma[i]);
    }

    __device__ __host__ bool check_aabb(const Ray &ray) const
    {
        float s[3]={ray.o.x, ray.o.y, ray.o.z}, v[3]={ray.d.x, ray.d.y, ray.d.z};
        Point3D pos;
        if (mi[0]<s[0]+eps && s[0]<ma[0]+eps && mi[1]<s[1]+eps && s[1]<ma[1]+eps && mi[2]<s[2]+eps && s[2]<ma[2]+eps)return true;
        for(int i=0;i<3;i++)
        {
            if (fabs(v[i]) > eps)
            {
                float t = ((v[i]>0?mi[i]:ma[i]) - s[i]) / v[i];
                if (t > eps)
                {
                    pos = ray.o + ray.d*t;
                    float p[3]={pos.x, pos.y, pos.z};
                    if (mi[(i+1)%3]<p[(i+1)%3]+eps && p[(i+1)%3]<ma[(i+1)%3]+eps && mi[(i+2)%3]<p[(i+2)%3]+eps && p[(i+2)%3]<ma[(i+2)%3]+eps)return true;
                }
            }
        }
        return false;
    }

};

struct KdTreeTemp
{
    size_t pos;
    bool in;
    float x;
    KdTreeTemp(size_t _pos, bool _in, float _x)
        :pos(_pos), in(_in), x(_x)
    {
    }

    bool operator < (const KdTreeTemp a) const
    {
        return x < a.x;
    }
};

struct KdTree
{
    KdTree *l, *r;
    KdTreeNode aabb;
    int w, cu_num;
    float splitl, splitr, split;
    std::vector<KdTreeNode> node;
    KdTreeNode *cu_node;

    bool check(const std::vector<Object> &objs, const Ray &ray, int &no, int &nv, float &dis);
    bool check_node(const std::vector<Object> &objs, const Ray &ray, int &no, int &nv, float &dis);
    KdTree(const std::vector<Object> &a, int s);
    KdTree(const std::vector<KdTreeNode> &a, const std::vector<KdTreeTemp> *com, int s);
    void create(const std::vector<KdTreeNode> &a, const std::vector<KdTreeTemp> *com, int s);

    ~KdTree()
    {
        if(!l)delete l;
        if(!r)delete r;
    }
};

#endif // KDTREE_H
