#ifndef KernelDefs
#define KernelDefs
#include <CL/cl.h>
#include "Kernels/KernelErrors.h"

enum KERNEL {
    UpdateZKernel,
    /*UpdateQKernel,*/
    /*GetNumFromPopsKernel,*/
    NumberOfKernels
};

enum BUFFER {
    QCL,
    PCL,
    ZCL,
    GENOCL,
    PREGENOCL,
    RECESSIVECL,
    NUMALLELESCL,
    POPFLAGCL,
    NUMAFROMPOPSCL,
    RANDCL, /* buffer for random numbers */
    ERRORCL, /* buffer for error codes */
    NumberOfBuffers
};
typedef struct CLDict {
    cl_kernel *kernels;
    cl_mem *buffers;
    cl_program program;
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_context context;
    cl_device_id device_id;
    cl_command_queue commands;
} CLDict;
#endif

