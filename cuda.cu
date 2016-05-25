#include <cstdio>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <curand_uniform.h>
#include <atomic>
#include <map>
#include "cuda_raytrace.h"
#include "polygon.h"
#include "object.h"
#include "sphere.h"
#include "kdtree.cuh"


#define CLEAR_MEM(p) if(p) {cudaFree(p); p=0;}

static Polygon *cu_polys;
static Sphere *cu_spheres;
static int numPoly=0, numSphere=0;
static unsigned char* cu_tex;
static unsigned char** cu_tex_pos;
static CuKdTree* cu_kdtrees;
static KdTreeNode* cu_nodes;
static curandState_t *rands;
static Point3D *ans;

__global__ void k_InitPoly(Polygon *points, int n, unsigned char* pos[], unsigned long long seed, curandState_t *rands)
{
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if(!idx)curand_init(seed, 0, 0, rands);
    if(idx >= n)return;
    Polygon &p = points[idx];
    const Point3D &a = p.points3d[0], &b=p.points3d[1], &c=p.points3d[2];
    p.nc = Point3D((b.y*c.z+c.y*a.z+a.y*b.z-b.y*a.z-a.y*c.z-c.y*b.z),
                (a.x*c.z+b.x*a.z+c.x*b.z-c.x*a.z-b.x*c.z-a.x*b.z),
                (a.x*b.y+b.x*c.y+c.x*a.y-c.x*b.y-b.x*a.y-a.x*c.y));
    p.nc.norm();
    p.tex = pos[idx];
}

int travKdTree(KdTree *tr, std::vector<CuKdTree> &trees, std::vector<KdTreeNode> &nodes, std::vector<int> &sum)
{
    if(!tr)return 0;
    CuKdTree tmp;
    tmp.aabb = tr->aabb;
    tmp.split = tr->split;
    tmp.w = tr->w;
    tmp.node = nodes.size();
    tmp.num = tr->node.size();
    for(int i=0;i<tr->node.size();i++)
    {
        KdTreeNode node = tr->node[i];
        node.nv = sum[node.no]+node.nv;
        nodes.push_back(node);
    }
    int t = trees.size();
    trees.push_back(tmp);
    int l = travKdTree(tr->l, trees, nodes, sum), r = travKdTree(tr->r, trees, nodes, sum);
    trees[t].l = l-t;
    trees[t].r = r-t;
    return t;
}

void copyToDevice(std::vector<Sphere> &spheres, std::vector<Object> &objs, KdTree *kdtree)
{
    for(int i=0;i<objs.size();i++)numPoly+=objs[i].polys.size();
    size_t size = numPoly * sizeof(Polygon);
    cudaMalloc((void**) &cu_polys, size);
    cudaMemset(cu_polys, 0, size);
    for(int i=0, sum=0;i<objs.size();sum+=objs[i].polys.size(), i++)
    {
        cudaMemcpy(cu_polys+sum, &objs[i].polys[0], objs[i].polys.size()*sizeof(Polygon), cudaMemcpyHostToDevice);
    }

    std::map<unsigned char*, unsigned char*> st;
    std::vector<unsigned char*> pos;
    int sumIm=0;
    for(int i=0;i<objs.size();i++)
    {
        for(int j=0;j<objs[i].polys.size();j++)
        {
            std::map<unsigned char*, unsigned char*>::iterator iter;
            if((iter = st.find(objs[i].polys[j].tex)) == st.end())
            {
                st[objs[i].polys[j].tex] = NULL;
                sumIm+=3*objs[i].polys[j].sizex*objs[i].polys[j].sizey;
            }
        }
    }
    st.clear();
    size = sumIm * sizeof(unsigned char);
    cudaMalloc((void**) &cu_tex, size);
    cudaMemset(cu_tex, 0, size);
    size = numPoly * sizeof(unsigned char*);
    cudaMalloc((void**) &cu_tex_pos, size);
    cudaMemset(cu_tex_pos, 0, size);
    unsigned char *cur = cu_tex;
    for(int i=0;i<objs.size();i++)
    {
        for(int j=0;j<objs[i].polys.size();j++)
        {
            std::map<unsigned char*, unsigned char*>::iterator iter;
            if((iter = st.find(objs[i].polys[j].tex)) == st.end())
            {
                st[objs[i].polys[j].tex] = cur;
                cudaMemcpy(cur, objs[i].polys[j].tex, 3*objs[i].polys[j].sizex*objs[i].polys[j].sizey*sizeof(unsigned char), cudaMemcpyHostToDevice);
                cur+=3*objs[i].polys[j].sizex*objs[i].polys[j].sizey;
                iter = st.find(objs[i].polys[j].tex);
            }
            pos.push_back(iter->second);
        }
    }
    cudaMemcpy(cu_tex_pos, &pos[0], numSphere*sizeof(unsigned char*), cudaMemcpyHostToDevice);

    dim3 grids((numSphere+31)/32);
    dim3 blocks(32);
    cudaMalloc((void**) &rands, sizeof(curandState_t));
    k_InitPoly<<<grids, blocks>>>(cu_polys, numSphere, cu_tex_pos, time(NULL), rands);
    cudaDeviceSynchronize();

    numSphere=spheres.size();
    size = numSphere * sizeof(Sphere);
    cudaMalloc((void**) &cu_spheres, size);
    cudaMemset(cu_spheres, 0, size);
    cudaMemcpy(cu_spheres, &spheres[0], size, cudaMemcpyHostToDevice);

    std::vector<CuKdTree> trees;
    std::vector<KdTreeNode> nodes;
    std::vector<int> sum;sum.push_back(0);
    for(int i=0;i<objs.size();i++)sum.push_back(sum.back()+objs[i].polys.size());
    travKdTree(kdtree, trees, nodes, sum);

    size = trees.size() * sizeof(CuKdTree);
    cudaMalloc((void**) &cu_kdtrees, size);
    cudaMemset(cu_kdtrees, 0, size);
    cudaMemcpy(cu_kdtrees, &trees[0], size, cudaMemcpyHostToDevice);

    size = nodes.size() * sizeof(KdTreeNode);
    cudaMalloc((void**) &cu_nodes, size);
    cudaMemset(cu_nodes, 0, size);
    cudaMemcpy(cu_nodes, &nodes[0], size, cudaMemcpyHostToDevice);

    cudaMalloc((void**) &ans, sizeof(Point3D));
    cudaMemset(ans, 0, sizeof(Point3D));

    cudaDeviceSynchronize();
}

__device__ bool cu_nodesCheck(const KdTreeNode *nodes, int n, const Polygon *polys, const Ray &ray, int &nv, double &dis)
{
    dis=1e20;
    for(int i=0;i<n;i++)
    {
        int x=nodes[i].nv;
        double t, u, v;
        if((t=polys[x].intersect(ray, u, v))<dis && t > eps)
            dis=t, nv=x;
    }
    return dis < 1e10;
}

__device__ bool cu_spheresCheck(const Sphere *spheres, int n, const Ray & r, int &id, double &t)
{
    id = n;
//    double d;
//    t = 1e20;
//    for (int i=0;i<n;i++)
//    {
//        if ((d = spheres[i].intersect(r)) && d < t)
//        {
//            t = d;
//            id = i;
//        }
//    }
    return t < 1e10;
}

__device__ bool cu_kdtreeCheck(const CuKdTree *tr, const KdTreeNode *nodes, const Polygon *polys, const Ray &ray, int &nv, double &dis)
{
    dis = 1e20;
    double dis1=1e20;
    if(!tr->aabb.check_aabb(ray))return false;
    double s[3]={ray.o.x, ray.o.y, ray.o.z}, v[3]={ray.d.x, ray.d.y, ray.d.z};
    int nv1, w=tr->w;
    bool t = false;
    if(s[w] < tr->split)
    {
        if(tr->l>0)t = cu_kdtreeCheck(tr+tr->l, nodes, polys, ray, nv, dis);
        cu_nodesCheck(nodes+tr->node, tr->num, polys, ray, nv1, dis1);
        if(dis1 < dis)nv=nv1, dis=dis1;
        if(tr->r>0 && !t && v[w]>0)cu_kdtreeCheck(tr+tr->r, nodes, polys, ray, nv1, dis1);
        if(dis1 < dis) nv=nv1, dis=dis1;
        return dis < 1e10;
    }
    if(s[w] > tr->split)
    {
        if(tr->r>0)t = cu_kdtreeCheck(tr+tr->r, nodes, polys, ray, nv, dis);
        cu_nodesCheck(nodes+tr->node, tr->num, polys, ray, nv1, dis1);
        if(dis1 < dis)nv=nv1, dis=dis1;
        if(tr->l>0 && !t && v[w]<0)cu_kdtreeCheck(tr+tr->l, nodes, polys, ray, nv1, dis1);
        if(dis1 < dis)nv=nv1, dis=dis1;
        return dis < 1e10;
    }
    cu_nodesCheck(nodes+tr->node, tr->num, polys, ray, nv, dis);
    if(v[w] < 0 && tr->l>0)cu_kdtreeCheck(tr+tr->r, nodes, polys, ray, nv1, dis1);
    if(v[w] > 0 && tr->r>0)cu_kdtreeCheck(tr+tr->l, nodes, polys, ray, nv1, dis1);
    if(dis1 < dis)nv=nv1, dis=dis1;
    return dis < 1e10;
}

__device__ Point3D cu_radiance(const CuKdTree *tr, const KdTreeNode *nodes, const Polygon *polys, const Sphere *spheres, const int &numSphere,
                               const Ray &r, int depth, bool into, curandState_t *rands)
{
    printf("111\n");
    double ts, to=1e20;
    int ids, idv;
    bool fs, fo;
    fo = cu_kdtreeCheck(tr, nodes, polys, r, idv, to);
    printf("222\n");
    fs = cu_spheresCheck(spheres, numSphere, r, ids, ts);
    printf("333\n");
    if(!fo && !fs)return Point3D();
    Point3D x, n, nl, f, lig;
    RType type;
//    if(fo)printf("%lf %lf %d %d\n", ts, to, ido, idv);
    if(ts < to || !fo)
    {
        const Sphere & obj = spheres[ids];
        x = r.o + r.d * ts, n = (x - obj.pos).norm(), nl = ((n*r.d) < 0 ? n : n * -1);
        f = obj.col, type = obj.type, lig = obj.lig;
    }
    else
    {
        double u, v;
        polys[idv].intersect(r, u, v);
        x = r.o + r.d * to;
        n = polys[idv].getn(x);
        nl = ((n*r.d) < 0 ? n : n * -1);
        f = polys[idv].getcol(u, v);
        type = polys[idv].type, lig = polys[idv].lig;
    }

    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;	// max refl
    if(++depth > 1)return lig;else if(depth > 5)f = f * (1 / p);

    if (type == DIFF)
    {
        double r1 = 2 * M_PI * curand_uniform(rands), r2 = curand_uniform(rands), r2s = sqrt(r2);
        Point3D w = nl, u = ((fabs(w.x) > .1 ? Point3D(0, 1, 0) : Point3D(1, 0, 0)) % w).norm(), v = w % u;
        Point3D d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return lig + f.mult(cu_radiance(tr, nodes, polys, spheres, numSphere, Ray(x-r.d*eps, d), depth, into, rands));
    } else if (type == SPEC)	// Ideal SPECULAR reflection
        return lig + f.mult(cu_radiance(tr, nodes, polys, spheres, numSphere, Ray(x-r.d*eps, r.d - n * 2 * (n*r.d)), depth, into, rands));

    Ray reflRay(x-r.d*eps, r.d - n * 2 * (n*r.d));	// Ideal dielectric REFRACTION
    double nc = 1, nt = 1.7, nnt = into ? nc / nt : nt / nc, ddn = r.d*nl, cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)	// Total internal reflection
        return lig + f.mult(cu_radiance(tr, nodes, polys, spheres, numSphere, reflRay, depth, into, rands));
    Point3D tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b),   c = 1 - (into ? -ddn : (tdir*n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return lig + f.mult(depth > 2 ? (curand_uniform(rands) < P ? cu_radiance(tr, nodes, polys, spheres, numSphere, reflRay, depth, into, rands) * RP :
                                                                 cu_radiance(tr, nodes, polys, spheres, numSphere, Ray(x+r.d*eps, tdir), depth, !into, rands) * TP) :
                                    cu_radiance(tr, nodes, polys, spheres, numSphere, reflRay, depth, into, rands) * Re
                                    + cu_radiance(tr, nodes, polys, spheres, numSphere, Ray(x+r.d*eps, tdir), depth, !into, rands) * Tr);
}

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
            (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

__global__ void k_radiance(const CuKdTree *tr, const KdTreeNode *nodes, const Polygon *polys, const Sphere *spheres, const int &numSphere,
                           int x, int y, int sizex, int sizey, int samps, Point3D cx, Point3D cy, Ray cam, curandState_t *rands, Point3D *ans)
{
    int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if(idx >= samps*4)return;
    int sx = idx&1, sy = idx&2;
    double r1 = 2 * curand_uniform(rands), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
    double r2 = 2 * curand_uniform(rands), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
    Point3D d = cx * (((sx + .5 + dx) / 2 + x) / sizex - .5) + cy * (((sy + .5 + dy) / 2 + y) / sizey - .5) + cam.d;
    Point3D t = cu_radiance(tr, nodes, polys, spheres, numSphere, Ray(cam.o + d * 140, d.norm()), 0, true, rands)*(.25 / samps);
    printf("%lf %lf %lf\n", t.x, t.y, t.z);
    atomicAdd(&(ans->x), t.x);
    atomicAdd(&(ans->y), t.y);
    atomicAdd(&(ans->z), t.z);
}

__global__ void test(const Sphere *spheres, const int &numSphere, Ray r, double *ans)
{
    int t;
    cu_spheresCheck(spheres, numSphere, r, t, *ans);
    printf("1231231%d\n", t);
}

Point3D cuda_radiance(int x, int y, int sizex, int sizey, int samps, Point3D cx, Point3D cy, Ray cam)
{
    dim3 grids((samps+7)/8);
    dim3 blocks(32);
    cudaMemset(ans, 0, sizeof(Point3D));
    k_radiance<<<grids, blocks>>>(cu_kdtrees, cu_nodes, cu_polys, cu_spheres, numSphere, x, y, sizex, sizey, samps, cx, cy, cam, rands, ans);
    cudaDeviceSynchronize();
    Point3D t;
    cudaMemcpy(&t, ans, sizeof(Point3D), cudaMemcpyDeviceToHost);
    return t;
}

bool cuda_test(std::vector<Sphere> &spheres, Ray r, double &dis)
{
//    cudaMemcpy(&spheres[0], cu_spheres, sizeof(Sphere)*spheres.size(), cudaMemcpyDeviceToHost);
//    cudaDeviceSynchronize();
    cudaMemset(ans, 0, sizeof(Point3D));
    test<<<dim3(), dim3()>>>(cu_spheres, numSphere, r, &ans->x);
    cudaDeviceSynchronize();
    Point3D t;
    cudaMemcpy(&t, ans, sizeof(Point3D), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    dis = t.x;
    return dis < 1e10;
}
