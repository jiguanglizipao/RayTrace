#ifndef KDTREE_CUH
#define KDTREE_CUH

#include <cuda_runtime.h>
#include "kdtree.h"
struct CuKdTree
{
    KdTreeNode aabb;
    int w, num;
    double split;
    int node, l, r;
};

#endif // KDTREE_H
