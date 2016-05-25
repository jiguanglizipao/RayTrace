#ifndef KDTREE_H
#define KDTREE_H

#include "point.h"
#include "polygon.h"
#include "object.h"
#include "algorithm"

struct KdTreeNode
{
    int no, nv;
    float mi[3], ma[3];
    KdTreeNode()
    {
        no=nv=-1;
    }

    KdTreeNode(const Polygon &a, int _no, int _nv)
        :no(_no), nv(_nv)
    {
        init();
        for(size_t i=0;i<3;i++)update(a.points3d[i]);
    }

    void init(float min=-1e10, float max=1e10)
    {
        mi[0]=mi[1]=mi[2]=max;
        ma[0]=ma[1]=ma[2]=min;
    }

    void update(const Point3D &a)
    {
        mi[0] = std::min(mi[0], a.x);
        mi[1] = std::min(mi[1], a.y);
        mi[2] = std::min(mi[2], a.z);
        ma[0] = std::max(ma[0], a.x);
        ma[1] = std::max(ma[1], a.y);
        ma[2] = std::max(ma[2], a.z);
    }

    void update(const KdTreeNode &a)
    {
        for(int i=0;i<3;i++)mi[i] = std::min(mi[i], a.mi[i]), ma[i] = std::max(ma[i], a.ma[i]);
    }

    bool check_aabb(Ray ray);

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
    int w;
    float splitl, splitr, split;
    std::vector<KdTreeNode> node;
    bool check(const std::vector<Object> &objs, Ray ray, int &no, int &nv, float &dis);
    bool check_node(const std::vector<Object> &objs, Ray ray, int &no, int &nv, float &dis);
    KdTree(const std::vector<Object> &a, int s);
    KdTree(const std::vector<KdTreeNode> &a, const std::vector<KdTreeTemp> *com, int s);
    void create(const std::vector<KdTreeNode> &a, const std::vector<KdTreeTemp> *com, int s);

    ~KdTree()
    {
        delete l;
        delete r;
    }
};

#endif // KDTREE_H
