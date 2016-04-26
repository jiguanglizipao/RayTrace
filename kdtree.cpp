#include "kdtree.h"
using namespace std;

KdTree::KdTree(const vector<Object> &a, int s)
    :l(NULL), r(NULL)
{
    vector<KdTreeNode>pre;
    for(size_t i=0;i<a.size();i++)for(size_t j=0;j<a[i].polys.size();j++)pre.push_back(KdTreeNode(a[i].polys[j], i, j));
    create(pre, s);
}

KdTree::KdTree(const vector<KdTreeNode> &a, int s)
    :l(NULL), r(NULL)
{
    create(a, s);
}

void KdTree::create(const vector<KdTreeNode> &pre, int s)
{
    w = s;
    if(pre.size() < 10)
    {
        node = pre;
        return;
    }
    double mi=1e10, ma=-1e10;
    for(size_t i=0;i<pre.size();i++)mi = min(mi, pre[i].mi[w]), ma = max(ma, pre[i].ma[w]);
    split = (ma+mi)*0.5;
    vector<KdTreeNode> vl, vr;
    vl.clear();vr.clear();node.clear();
    for(size_t i=0;i<pre.size();i++)
    {
        if(pre[i].ma[w] < split)
            vl.push_back(pre[i]);
        else if(pre[i].mi[w] > split)
            vr.push_back(pre[i]);
        else
            node.push_back(pre[i]);
    }
    if(!vl.empty())l = new KdTree(vl, (w+1)%3);
    if(!vr.empty())r = new KdTree(vr, (w+1)%3);
}

bool KdTree::check_node(const std::vector<Object> &objs, Point3D view, Point3D ray, int &no, int &nv, Point3D &p, double &dis)
{
    dis=1.1e10;
    double raydis = sqrt(ray*ray);
    for(size_t l=0;l<node.size();l++)
    {
        int i=node[l].no, j=node[l].nv;
        double k = -(objs[i].polys[j].d+objs[i].polys[j].n*view)/(objs[i].polys[j].n*ray);
        if(k < eps)continue;
        Point3D tmp = view+k*ray;
        if(objs[i].polys[j].checkInside(tmp) == Polygon::inside)
            if(k*raydis < dis)dis=k*raydis, no=i, nv=j, p=tmp;
    }
    return dis < 1e10;
}

bool KdTree::check(const std::vector<Object> &objs, Point3D view, Point3D ray, int &no, int &nv, Point3D &p, double &dis)
{
    double s[3]={view.x, view.y, view.z}, v[3]={ray.x, ray.y, ray.z};
    int no1, nv1;
    double dis1=1.1e10;dis = 1.1e10;
    Point3D p1;
    if(s[w] < split)
    {
        if(l)l->check(objs, view, ray, no, nv, p, dis);
        check_node(objs, view, ray, no1, nv1, p1, dis1);
        if(dis1 < dis)no=no1, nv=nv1, p=p1, dis=dis1;
        if(dis < 1e10)return true;else if(v[w] < 0)return false;
        if(r)r->check(objs, view, ray, no, nv, p, dis);
        return dis < 1e10;
    }
    if(s[w] > split)
    {
        if(r)r->check(objs, view, ray, no, nv, p, dis);
        check_node(objs, view, ray, no1, nv1, p1, dis1);
        if(dis1 < dis)no=no1, nv=nv1, p=p1, dis=dis1;
        if(dis < 1e10)return true;else if(v[w] > 0)return false;
        if(l)l->check(objs, view, ray, no, nv, p, dis);
        return dis < 1e10;
    }
    check_node(objs, view, ray, no1, nv1, p1, dis1);
    if(v[w] < 0 && l)l->check(objs, view, ray, no, nv, p, dis);
    if(v[w] > 0 && r)r->check(objs, view, ray, no, nv, p, dis);
    if(dis1 < dis)no=no1, nv=nv1, p=p1, dis=dis1;
    return dis < 1e10;
}
