#ifndef CUDA_H
#define CUDA_H
#include <cuda_runtime.h>
#include <cuda.h>
#include <vector>
#include "polygon.h"
#include "object.h"
#include "sphere.h"
#include "kdtree.cuh"
int check_device();
void print_device_info(int rank, int dev);
void copyToDevice(std::vector<Sphere> &spheres, std::vector<Object> &objs, KdTree *kdtree);
Point3D cuda_radiance(int x, int y, int sizex, int sizey, int samps, Point3D cx, Point3D cy, Ray cam);
#endif
