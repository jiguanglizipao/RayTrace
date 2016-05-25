#ifndef CUDA_H
#define CUDA_H
#include <cuda_runtime.h>
#include <cuda.h>
#include "polygon.h"
const int GPU_PER_NODE = 1;
int check_device();
void print_device_info(int rank, int dev);
#endif
