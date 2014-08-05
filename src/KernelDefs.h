#ifndef KernelDefs
#define KernelDefs
#include <CL/cl.h>
#include "Kernels/KernelErrors.h"
/*TODO: Determine this number, used for reductions */
#define MAXGROUPS 31

enum KERNEL {
    UpdateZKernel,
    GetNumFromPopsKernel,
    UpdatePKernel,
    mapReduceLogDiffsKernel,
    RDirichletSampleKernel,
    MetroAcceptTestKernel,
    GetNumLociPopsKernel,
    UpdQDirichletKernel,
    FillArrayWRandomKernel,
    InitRandGenKernel,
    UpdateFstKernel,
    PopNormals,
    UpdateAlphaKernel,
    NonIndUpdateEpsilonKernel,
    DataCollectPopKernel,
    DataCollectIndKernel,
    DataCollectLocKernel,
    CalcLikeKernel,
    ComputeProbFinishKernel,
    MapReduceLogLikeKernel,
    NumberOfKernels
};


enum BUFFER {
    QCL,
    QSUMCL,
    PCL,
    PSUMCL,
    ZCL,
    GENOCL,
    PREGENOCL,
    RECESSIVECL,
    NUMALLELESCL,
    LAMBDACL,
    LAMBDASUMCL,
    POPFLAGCL,
    NUMAFROMPOPSCL,
    NUMLOCIPOPSCL,
    EPSILONCL,
    EPSILONSUMCL,
    FSTCL,
    FSTSUMCL,
    NORMSCL,
    /*LOGTERMSCL,*/
    RANDGENSCL,
    LOGDIFFSCL,
    REDUCERESULTSCL,
    TESTQCL,
    ALPHACL,
    ALPHASUMCL,
    ANCESTDISTCL,
    RANDCL, /* buffer for random numbers */
    ERRORCL, /* buffer for error codes */
    LIKECL,
    SUMLIKESCL,
    SUMSQLIKESCL,
    LOGLIKESCL,
    INDLIKECL,
    SUMINDLIKECL,
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

