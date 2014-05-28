#include <CL/cl.h>
typedef struct CLDict {
    int numProgramsInDict;
    char **ProgramNames;
    cl_program *Programs;    
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_context context;
    cl_device_id device_id; 
    cl_command_queue commands;
} CLDict;

extern int InitCLDict(CLDict *clDictToInit);
extern void ReleaseCLDict(CLDict *clDict);
extern int CompileProgram(char * programFilename, CLDict *clDict);
