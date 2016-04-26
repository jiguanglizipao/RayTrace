#ifndef KDTREE_H
#define KDTREE_H

#include "point.h"
#include "polygon.h"
#include "object.h"
#include "algorithm"

struct KdTreeNode
{
    int no, nv;
    double mi[3], ma[3];
    KdTreeNode(const Polygon &a, int _no, int _nv)
        :no(_no), nv(_nv)
    {
        mi[0]=mi[1]=mi[2]=1e10;
        ma[0]=ma[1]=ma[2]=-1e10;
        for(size_t i=0;i<a.points3d.size();i++)
        {
            mi[0] = std::min(mi[0], a.points3d[i].x);
            mi[1] = std::min(mi[1], a.points3d[i].y);
            mi[2] = std::min(mi[2], a.points3d[i].z);
            ma[0] = std::max(ma[0], a.points3d[i].x);
            ma[1] = std::max(ma[1], a.points3d[i].y);
            ma[2] = std::max(ma[2], a.points3d[i].z);
        }
    }
};

struct KdTree
{
    KdTree *l, *r;
    int w;
    double splitl, splitr, split;
    std::vector<KdTreeNode> node;
    bool check(const std::vector<Object> &objs, Point3D view, Point3D ray, int &no, int &nv, Point3D &p, double &dis);
    bool check_node(const std::vector<Object> &objs, Point3D view, Point3D ray, int &no, int &nv, Point3D &p, double &dis);
    KdTree(const std::vector<Object> &a, int s);
    KdTree(const std::vector<KdTreeNode> &a, int s);
    void create(const std::vector<KdTreeNode> &a, int s);
    ~KdTree()
    {
        delete l;
        delete r;
    }
};

#endif // KDTREE_H
