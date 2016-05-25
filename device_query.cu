#include "cuda.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
int check_device()
{
    int deviceCount;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if(err != cudaSuccess)
    {
      printf("cudaGetDeviceCount returned: %s\n", cudaGetErrorString(err));
    }
    return deviceCount;
}

void print_device_info(int rank, int dev)
{
   cudaDeviceProp deviceProp;
   cudaGetDeviceProperties(&deviceProp, dev);
  
    printf("\nRank %d - Device %d: \"%s\"\n", rank, dev, deviceProp.name);
    printf("  Major revision number:                         %d\n",
               deviceProp.major);
    printf("  Minor revision number:                         %d\n",
               deviceProp.minor);
    printf("  Total amount of global memory:                 %lu bytes\n",
               deviceProp.totalGlobalMem);
#if CUDART_VERSION >= 2000
    printf("  Number of multiprocessors:                     %d\n",
               deviceProp.multiProcessorCount);
    printf("  Number of cores:                               %d\n",
               8 * deviceProp.multiProcessorCount);
#endif
    printf("  Total amount of constant memory:               %lu bytes\n",
               deviceProp.totalConstMem); 
    printf("  Total amount of shared memory per block:       %lu bytes\n",
               deviceProp.sharedMemPerBlock);
    printf("  Total number of registers available per block: %d\n",
               deviceProp.regsPerBlock);
    printf("  Warp size:                                     %d\n",
               deviceProp.warpSize);
    printf("  Maximum number of threads per block:           %d\n",
               deviceProp.maxThreadsPerBlock);
    printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
    printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
    printf("  Maximum memory pitch:                          %lu bytes\n",
               deviceProp.memPitch);
    printf("  Texture alignment:                             %lu bytes\n",
               deviceProp.textureAlignment);
    printf("  Clock rate:                                    %.2f GHz\n",
               deviceProp.clockRate * 1e-6f);
#if CUDART_VERSION >= 2000
    printf("  Concurrent copy and execution:                 %s\n",
               deviceProp.deviceOverlap ? "Yes" : "No");
#endif
    if(!deviceProp.deviceOverlap)
	printf("Device will not handle overlaps, so no speed up from streams \n");
}
