#ifndef KernelDefs
#define KernelDefs
#include <CL/cl.h>
#include "Kernels/KernelErrors.h"

enum KERNEL {
    UpdateZKernel,
    GetNumFromPopsKernel,
    UpdatePKernel,
    mapLogDiffsKernel,
    reduceLogDiffsKernel,
    NumberOfKernels
};


enum BUFFER {
    QCL,
    PCL,
    LOGPCL,
    ZCL,
    GENOCL,
    PREGENOCL,
    RECESSIVECL,
    NUMALLELESCL,
    LAMBDACL,
    POPFLAGCL,
    NUMAFROMPOPSCL,
    EPSILONCL,
    FSTCL,
    LOGTERMSCL,
    LOGDIFFSCL,
    TESTQCL,
    RANDCL, /* buffer for random numbers */
    ERRORCL, /* buffer for error codes */
    NumberOfBuffers
};
typedef struct CLDict {
    cl_kernel *kernels;
    cl_mem *buffers;
    cl_program program;
    size_t *locals;
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_context context;
    cl_device_id device_id;
    cl_command_queue commands;
} CLDict;
#endif

