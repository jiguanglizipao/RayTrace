#include "kdtree.h"
#include <cstdio>
#include <algorithm>
using namespace std;

KdTree::KdTree(const vector<Object> &a, int s)
    :l(NULL), r(NULL)
{
    vector<KdTreeNode>pre;
    vector<KdTreeTemp>com[3];
    for(int i=0;i<3;i++)com[i].clear();
    for(size_t i=0;i<a.size();i++)for(size_t j=0;j<a[i].polys.size();j++)
    {
        pre.push_back(KdTreeNode(a[i].polys[j], i, j));
        for(int k=0;k<3;k++)
        {
            com[k].push_back(KdTreeTemp(pre.size()-1, true, pre.back().mi[k]));
            com[k].push_back(KdTreeTemp(pre.size()-1, false, pre.back().ma[k]));
        }
    }
    for(int i=0;i<3;i++)sort(com[i].begin(), com[i].end());
    create(pre, com, s);
}

KdTree::KdTree(const vector<KdTreeNode> &a, const std::vector<KdTreeTemp> *com, int s)
    :l(NULL), r(NULL)
{
    create(a, com, s);
}

void KdTree::create(const vector<KdTreeNode> &pre, const std::vector<KdTreeTemp> *com, int s)
{
    w = s;
    aabb.init();
    for(size_t i=0;i<pre.size();i++)aabb.update(pre[i]);
    if(pre.size() < 16)
    {
        node = pre;
        return;
    }

//    int suml=0, summ=0, sumr=pre.size();
//    double mi = 0.5*pre.size()+pow(max(suml, sumr), 2.0/3.0);
//    split = com[w][0].x-eps;
//    for(size_t i=0;i<com[w].size();)
//    {
//        double c = com[w][i].x;
//        for(;i<com[w].size() && com[w][i].x<c+eps;i++)
//        {
//            if(com[w][i].in)summ++, sumr--;else suml++, summ--;
//        }
//        double t = 0.5*abs(suml-sumr)+pow(max(suml, sumr), 2.0/3.0)+summ;
//        if(t < mi)mi = t, split = c;
//    }

//    vector<KdTreeNode> vl, vr;
//    vector<KdTreeTemp> coml[3], comr[3];
//    vector<pair<int, int> > flag;
//    for(size_t i=0;i<pre.size();i++)
//    {
//        if(pre[i].ma[w] < split)
//            vl.push_back(pre[i]), flag.push_back(make_pair(0, vl.size()-1));
//        else if(pre[i].mi[w] > split)
//            vr.push_back(pre[i]), flag.push_back(make_pair(1, vr.size()-1));
//        else
//            node.push_back(pre[i]), flag.push_back(make_pair(2, node.size()-1));
//    }
//    for(int i=0;i<3;i++)
//    {
//        coml[i].clear();comr[i].clear();
//        for(size_t j=0;j<com[i].size();j++)
//        {
//            if(flag[com[i][j].pos].first == 0)
//                coml[i].push_back(KdTreeTemp(flag[com[i][j].pos].second, com[i][j].in, com[i][j].x));
//            if(flag[com[i][j].pos].first == 1)
//                comr[i].push_back(KdTreeTemp(flag[com[i][j].pos].second, com[i][j].in, com[i][j].x));
//        }
//    }
//    if(!vl.empty())l = new KdTree(vl, coml, (w+1)%3);
//    if(!vr.empty())r = new KdTree(vr, comr, (w+1)%3);

    split = (aabb.mi[w]+aabb.ma[w])*0.5;
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
    if(!vl.empty())l = new KdTree(vl, com, (w+1)%3);
    if(!vr.empty())r = new KdTree(vr, com, (w+1)%3);
}

bool KdTree::check_node(const std::vector<Object> &objs, const Ray &ray, int &no, int &nv, double &dis)
{
    dis=1e20;
    for(size_t l=0;l<node.size();l++)
    {
        int i=node[l].no, j=node[l].nv;
        double t, u, v;
        if((t=objs[i].polys[j].intersect(ray, u, v))<dis && t > eps)
            dis=t, no=i, nv=j;
    }
    return dis < 1e10;
}

bool KdTree::check(const std::vector<Object> &objs, const Ray &ray, int &no, int &nv, double &dis)
{
    dis = 1e20;
    double dis1=1e20;
    if(!aabb.check_aabb(ray))return false;
    double s[3]={ray.o.x, ray.o.y, ray.o.z}, v[3]={ray.d.x, ray.d.y, ray.d.z};
    int no1, nv1;
    bool t = false;
    if(s[w] < split)
    {
        if(l)t = l->check(objs, ray, no, nv, dis);
        check_node(objs, ray, no1, nv1, dis1);
        if(dis1 < dis)no=no1, nv=nv1, dis=dis1;
        if(r && !t && v[w]>0)r->check(objs, ray, no1, nv1, dis1);
        if(dis1 < dis)no=no1, nv=nv1, dis=dis1;
        return dis < 1e10;
    }
    if(s[w] > split)
    {
        if(r)t = r->check(objs, ray, no, nv, dis);
        check_node(objs, ray, no1, nv1, dis1);
        if(dis1 < dis)no=no1, nv=nv1, dis=dis1;
        if(l && !t && v[w]<0)l->check(objs, ray, no1, nv1, dis1);
        if(dis1 < dis)no=no1, nv=nv1, dis=dis1;
        return dis < 1e10;
    }
    check_node(objs, ray, no1, nv1, dis1);
    if(v[w] < 0 && l)l->check(objs, ray, no, nv, dis);
    if(v[w] > 0 && r)r->check(objs, ray, no, nv, dis);
    if(dis1 < dis)no=no1, nv=nv1, dis=dis1;
    return dis < 1e10;
}
