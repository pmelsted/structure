#include <CL/cl.h>

enum KERNEL {
    UpdateZKernel,
    /*UpdateQKernel,*/
    NumberOfKernels
};

typedef struct CLDict {
    int numKernelsInDict;
    cl_kernel *kernels;
    cl_program program;    
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_context context;
    cl_device_id device_id; 
    cl_command_queue commands;
} CLDict;
