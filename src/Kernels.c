#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>
#include "KernelDefs.h"
#include "structure.h"

#define MAX_SOURCE_SIZE (0x100000)
#define USEGPU 1


void printCLErr(cl_int err)
{
    switch (err) {
    case CL_SUCCESS:
        printf("CL_SUCCES S\n");
        break;
    case CL_DEVICE_NOT_FOUND:
        printf("CL_DEVICE_NOT_FOUND\n");
        break;
    case CL_DEVICE_NOT_AVAILABLE:
        printf("CL_DEVICE_NOT_AVAILABLE\n");
        break;
    case CL_COMPILER_NOT_AVAILABLE:
        printf("CL_COMPILER_NOT_AVAILABLE\n");
        break;
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:
        printf("CL_MEM_OBJECT_ALLOCATION_FAILURE\n");
        break;
    case CL_OUT_OF_RESOURCES:
        printf("CL_OUT_OF_RESOURCES\n");
        break;
    case CL_OUT_OF_HOST_MEMORY:
        printf("CL_OUT_OF_HOST_MEMORY\n");
        break;
    case CL_PROFILING_INFO_NOT_AVAILABLE:
        printf("CL_PROFILING_INFO_NOT_AVAILABLE\n");
        break;
    case CL_MEM_COPY_OVERLAP:
        printf("CL_MEM_COPY_OVERLAP\n");
        break;
    case CL_IMAGE_FORMAT_MISMATCH:
        printf("CL_IMAGE_FORMAT_MISMATCH\n");
        break;
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:
        printf("CL_IMAGE_FORMAT_NOT_SUPPORTED\n");
        break;
    case CL_BUILD_PROGRAM_FAILURE:
        printf("CL_BUILD_PROGRAM_FAILURE\n");
        break;
    case CL_MAP_FAILURE:
        printf("CL_MAP_FAILURE\n");
        break;
    case CL_INVALID_VALUE:
        printf("CL_INVALID_VALUE\n");
        break;
    case CL_INVALID_DEVICE_TYPE:
        printf("CL_INVALID_DEVICE_TYPE\n");
        break;
    case CL_INVALID_PLATFORM:
        printf("CL_INVALID_PLATFORM\n");
        break;
    case CL_INVALID_DEVICE:
        printf("CL_INVALID_DEVICE\n");
        break;
    case CL_INVALID_CONTEXT:
        printf("CL_INVALID_CONTEXT\n");
        break;
    case CL_INVALID_QUEUE_PROPERTIES:
        printf("CL_INVALID_QUEUE_PROPERTIES\n");
        break;
    case CL_INVALID_COMMAND_QUEUE:
        printf("CL_INVALID_COMMAND_QUEUE\n");
        break;
    case CL_INVALID_HOST_PTR:
        printf("CL_INVALID_HOST_PTR\n");
        break;
    case CL_INVALID_MEM_OBJECT:
        printf("CL_INVALID_MEM_OBJECT\n");
        break;
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
        printf("CL_INVALID_IMAGE_FORMAT_DESCRIPTOR\n");
        break;
    case CL_INVALID_IMAGE_SIZE:
        printf("CL_INVALID_IMAGE_SIZE\n");
        break;
    case CL_INVALID_SAMPLER:
        printf("CL_INVALID_SAMPLER\n");
        break;
    case CL_INVALID_BINARY:
        printf("CL_INVALID_BINARY\n");
        break;
    case CL_INVALID_BUILD_OPTIONS:
        printf("CL_INVALID_BUILD_OPTIONS\n");
        break;
    case CL_INVALID_PROGRAM:
        printf("CL_INVALID_PROGRAM\n");
        break;
    case CL_INVALID_PROGRAM_EXECUTABLE:
        printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
        break;
    case CL_INVALID_KERNEL_NAME:
        printf("CL_INVALID_KERNEL_NAME\n");
        break;
    case CL_INVALID_KERNEL_DEFINITION:
        printf("CL_INVALID_KERNEL_DEFINITION\n");
        break;
    case CL_INVALID_KERNEL:
        printf("CL_INVALID_KERNEL\n");
        break;
    case CL_INVALID_ARG_INDEX:
        printf("CL_INVALID_ARG_INDEX\n");
        break;
    case CL_INVALID_ARG_VALUE:
        printf("CL_INVALID_ARG_VALUE\n");
        break;
    case CL_INVALID_ARG_SIZE:
        printf("CL_INVALID_ARG_SIZE\n");
        break;
    case CL_INVALID_KERNEL_ARGS:
        printf("CL_INVALID_KERNEL_ARGS\n");
        break;
    case CL_INVALID_WORK_DIMENSION:
        printf("CL_INVALID_WORK_DIMENSION\n");
        break;
    case CL_INVALID_WORK_GROUP_SIZE:
        printf("CL_INVALID_WORK_GROUP_SIZE\n");
        break;
    case CL_INVALID_WORK_ITEM_SIZE:
        printf("CL_INVALID_WORK_ITEM_SIZE\n");
        break;
    case CL_INVALID_GLOBAL_OFFSET:
        printf("CL_INVALID_GLOBAL_OFFSET\n");
        break;
    case CL_INVALID_EVENT_WAIT_LIST:
        printf("CL_INVALID_EVENT_WAIT_LIST\n");
        break;
    case CL_INVALID_EVENT:
        printf("CL_INVALID_EVENT\n");
        break;
    case CL_INVALID_OPERATION:
        printf("CL_INVALID_OPERATION\n");
        break;
    case CL_INVALID_GL_OBJECT:
        printf("CL_INVALID_GL_OBJECT\n");
        break;
    case CL_INVALID_BUFFER_SIZE:
        printf("CL_INVALID_BUFFER_SIZE\n");
        break;
    case CL_INVALID_MIP_LEVEL:
        printf("CL_INVALID_MIP_LEVEL\n");
        break;
    case CL_INVALID_GLOBAL_WORK_SIZE:
        printf("CL_INVALID_GLOBAL_WORK_SIZE\n");
        break;
    }
}

void PrintKernelError(int error)
{
    switch(error) {
    case KERNEL_SUCCESS:
        printf("No error, Kernel success\n");
        break;
    case KERNEL_OUT_OF_BOUNDS:
        printf("Kernel out of bounds\n");
        break;
    default:
        printf("Unknown error %d\n",error);
    }
}

void ReleaseCLDict(CLDict *clDict)
{
    int i;
    clReleaseProgram(clDict->program);
    for(i = 0; i < NumberOfKernels; ++i) {
        clReleaseKernel(clDict->kernels[i]);
    }
    free(clDict->kernels);
    for(i = 0; i < NumberOfBuffers; ++i) {
        clReleaseMemObject(clDict->buffers[i]);
    }
    free(clDict->buffers);
    free(clDict->locals);
    clReleaseCommandQueue(clDict->commands);
    clReleaseContext(clDict->context);
    free(clDict);
}

void handleCLErr(cl_int err,CLDict *clDict, char * message)
{
    if (err != CL_SUCCESS) {
        printf("CL Error:\n");
        printf("%s\n",message);
        printCLErr(err);
        ReleaseCLDict(clDict);
        exit(EXIT_FAILURE);
    }
}

void setKernelArg(CLDict *clDict, enum KERNEL kernel, enum BUFFER buffer,int argnum){
    cl_int err;
    err = 0;
    /*=========================== Update Z =================*/
    err  = clSetKernelArg(clDict->kernels[kernel], argnum, sizeof(cl_mem),
                          &(clDict->buffers[buffer]));
    handleCLErr(err, clDict,"Failed to set arg!");

}


void setKernelArgNULL(CLDict *clDict, enum KERNEL kernel,size_t size, void *arg, int argnum){
    cl_int err;
    err = 0;
    err  = clSetKernelArg(clDict->kernels[kernel], argnum, size, arg);
    handleCLErr(err, clDict,"Failed to set arg!");

}

void setKernelArgs(CLDict *clDict)
{
    cl_int err;
    err = 0;
    /*=========================== Update Z =================*/
    setKernelArg(clDict,UpdateZKernel,QCL,0);
    setKernelArg(clDict,UpdateZKernel,PCL,1);
    setKernelArg(clDict,UpdateZKernel,GENOCL,2);
    setKernelArg(clDict,UpdateZKernel,RANDCL,3);
    setKernelArg(clDict,UpdateZKernel,ZCL,4);
    setKernelArg(clDict,UpdateZKernel,ERRORCL,5);
    /*=====================================================*/



    /* ========== Get num from pops ======= */
    setKernelArg(clDict,GetNumFromPopsKernel,ZCL,0);
    setKernelArg(clDict,GetNumFromPopsKernel,GENOCL,1);
    setKernelArg(clDict,GetNumFromPopsKernel,NUMALLELESCL,2);
    setKernelArg(clDict,GetNumFromPopsKernel,POPFLAGCL,3);
    setKernelArg(clDict,GetNumFromPopsKernel,NUMAFROMPOPSCL,4);
    setKernelArg(clDict,GetNumFromPopsKernel,ERRORCL,5);

    /*=====================================================*/


    /*=========================== Update P =================*/
    setKernelArg(clDict,UpdatePKernel,PCL,0);
    setKernelArg(clDict,UpdatePKernel,LOGPCL,1);
    setKernelArg(clDict,UpdatePKernel,NUMALLELESCL,2);
    setKernelArg(clDict,UpdatePKernel,NUMAFROMPOPSCL,3);
    setKernelArg(clDict,UpdatePKernel,RANDCL,4);
    setKernelArg(clDict,UpdatePKernel,ERRORCL,5);


    if (FREQSCORR) {
        setKernelArg(clDict,UpdatePKernel,EPSILONCL,6);
        setKernelArg(clDict,UpdatePKernel,FSTCL,7);
    } else {
        setKernelArg(clDict,UpdatePKernel,LAMBDACL,6);
    }
    /*=====================================================*/


    /*===== map log diffs *==== */

    /*setKernelArg(clDict,mapLogDiffsKernel,LOGTERMSCL,0);*/
    /*setKernelArg(clDict,mapLogDiffsKernel,QCL,1);*/
    /*setKernelArg(clDict,mapLogDiffsKernel,TESTQCL,2);*/
    /*setKernelArg(clDict,mapLogDiffsKernel,PCL,3);*/
    /*setKernelArg(clDict,mapLogDiffsKernel,GENOCL,4);*/
    /*setKernelArg(clDict,mapLogDiffsKernel,ERRORCL,5);*/

    /*=====  reduce log diffs *==== */
    /*setKernelArg(clDict,reduceLogDiffsKernel,LOGTERMSCL,0);*/
    /*setKernelArg(clDict,reduceLogDiffsKernel,LOGDIFFSCL,1);*/
    /*setKernelArgNULL(clDict,reduceLogDiffsKernel,sizeof(double)*NUMLOCI,NULL,2);*/

    /*=====  mapreduce log diffs *==== */
    setKernelArg(clDict,mapReduceLogDiffsKernel,QCL,0);
    setKernelArg(clDict,mapReduceLogDiffsKernel,TESTQCL,1);
    setKernelArg(clDict,mapReduceLogDiffsKernel,PCL,2);
    setKernelArg(clDict,mapReduceLogDiffsKernel,GENOCL,3);
    setKernelArg(clDict,mapReduceLogDiffsKernel,LOGDIFFSCL,4);
    setKernelArgNULL(clDict,mapReduceLogDiffsKernel,sizeof(double)*NUMLOCI,NULL,5);
    /*setKernelArgNULL(clDict,mapReduceLogDiffsKernel,sizeof(double)*NUMLOCI,NULL,5);*/

    /* RDirichlet sample */
    setKernelArg(clDict,RDirichletSampleKernel,ALPHACL,0);
    setKernelArg(clDict,RDirichletSampleKernel,RANDCL,1);
    setKernelArg(clDict,RDirichletSampleKernel,TESTQCL,2);

    /* Metro accept test*/
    setKernelArg(clDict,MetroAcceptTestKernel,TESTQCL,0);
    setKernelArg(clDict,MetroAcceptTestKernel,QCL,1);
    setKernelArg(clDict,MetroAcceptTestKernel,RANDCL,2);
    setKernelArg(clDict,MetroAcceptTestKernel,LOGDIFFSCL,3);
    setKernelArg(clDict,MetroAcceptTestKernel,POPFLAGCL,4);

    /* Get num loci pops */
    setKernelArg(clDict,GetNumLociPopsKernel,ZCL,0);
    setKernelArg(clDict,GetNumLociPopsKernel,POPFLAGCL,1);
    setKernelArg(clDict,GetNumLociPopsKernel,NUMLOCIPOPSCL,2);

    /* UpdQ dirichlet */
    setKernelArg(clDict,UpdQDirichletKernel,ALPHACL,0);
    setKernelArg(clDict,UpdQDirichletKernel,NUMLOCIPOPSCL,1);
    setKernelArg(clDict,UpdQDirichletKernel,RANDCL,2);
    setKernelArg(clDict,UpdQDirichletKernel,QCL,3);

}

void getLocal(CLDict *clDict,enum KERNEL kernel){
    cl_int err;
    err = clGetKernelWorkGroupInfo(clDict->kernels[kernel],
                                   clDict->device_id,
                                   CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
                                   &(clDict->locals[kernel]), NULL);
    handleCLErr(err, clDict, "Error: Failed to retrieve kernel work group info!");
}


void getLocals(CLDict *clDict)
{
    getLocal(clDict,UpdateZKernel);
    getLocal(clDict,GetNumFromPopsKernel);
    getLocal(clDict,UpdatePKernel);

}

char * searchReplace(char * string,  char *toReplace[], char *replacements[],
                     int numReplacements)
{
    int i = 0;
    char *locOfToRep;
    char *toRep;
    char *rep;
    int lenToRep,lenStr,lenAfterLocRep;
    static char buffer[MAX_SOURCE_SIZE];
    strncpy(buffer,string,MAX_SOURCE_SIZE);
    for(i = 0; i < numReplacements; ++i) {
        toRep = toReplace[i];
        rep = replacements[i];
        /*if str not in the string, skip it */
        if (!(locOfToRep = strstr(buffer,toRep))) {
            printf("Warning: %s not found in string!\n",toRep);
            continue;
        }
        lenToRep = strlen(toRep);
        lenStr = strlen(buffer);
        lenAfterLocRep = strlen(locOfToRep);
        /*Print the string upto the pointer, then the val, and then the rest of the string.*/
        sprintf(buffer, "%.*s%s%s", lenStr-lenAfterLocRep, buffer,rep,
                locOfToRep+lenToRep);
        /*Will break if buffer is longer than string, so we restrict ourselves to the longest length.*/
    }
    return buffer;
}

void preProcessSource(char * kernelSource, size_t *source_size,
                      char *names[], char *vals[], int numVals)
{
    char * processedSource;
    processedSource = searchReplace(kernelSource, names, vals, numVals);
    strncpy(kernelSource,processedSource,MAX_SOURCE_SIZE);
    *source_size = strlen(kernelSource);
}



int initKernel(CLDict *clDict,char * kernelName, enum KERNEL kernelEnumVal)
{
    cl_kernel kernel;
    int err;
    printf("%s, %d\n",kernelName,kernelEnumVal);
    kernel = clCreateKernel(clDict->program, kernelName, &err);
    if (!kernel || err != CL_SUCCESS) {
        switch(err) {
        case CL_INVALID_PROGRAM:
            printf("invalid program\n");
            break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
            printf("invalid exec\n");
            break;
        case CL_INVALID_KERNEL_NAME:
            printf("invalid name\n");
            break;
        case CL_INVALID_KERNEL_DEFINITION:
            printf("invalid def\n");
            break;
        case CL_INVALID_VALUE:
            printf("invalid val\n");
            break;
        case CL_OUT_OF_HOST_MEMORY:
            printf("invalid mem\n");
            break;
        case CL_SUCCESS:
            printf("no kernel\n");
            break;
        default:
            printf("%d\n", err);
        }
        printf("Error: Failed to create compute kernel %s!\n", kernelName);
        return EXIT_FAILURE;
    }
    clDict->kernels[kernelEnumVal] = kernel;
    return EXIT_SUCCESS;
}


/*
 * compiles the program with filename programFilename, and replaces the names in names with the values in vals.
 */
int CompileKernels(CLDict *clDict,  char *options)
{
    FILE *fp;
    char *KernelSource;
    size_t source_size;

    int err;
    int i;
    cl_program program;
    /*cl_uint numkernels;*/

    /*cl_int ret;*/


    char *KERNELNAMES[NumberOfKernels] = {"UpdateZ","GetNumFromPops","UpdateP","mapReduceLogDiffs","Dirichlet", "MetroAcceptTest","GetNumLociPops","UpdQDirichlet"};

    /* Load the source code containing the kernels*/
    fp = fopen("Kernels/Kernels.cl", "r");
    if (!fp) {
        fprintf(stderr, "Failed to load kernel file Kernels.cl\n");
        return EXIT_FAILURE;
    }


    KernelSource = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(KernelSource, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    /*preProcessSource(KernelSource, &source_size, names,vals,numVals);*/
    program = clCreateProgramWithSource(clDict->context, 1,
                                        (const char **) & KernelSource, NULL, &err);
    if (!program) {
        printf("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }


    free(KernelSource);


    err = clBuildProgram(program, 1, &clDict->device_id, options, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t len;
        char buffer[32768];
        printCLErr(err);
        printf("Error: Failed to build program executable!\n");
        printf("Options:");
        printf("%s\n",options);
        printf("Error:");
        printCLErr(err);
        clGetProgramBuildInfo(program, clDict->device_id, CL_PROGRAM_BUILD_LOG,
                              sizeof(buffer),buffer, &len);
        printf("%s\n", buffer);
        printf("%d\n", (int) len);
        return EXIT_FAILURE;
    }

    clDict->program = program;
    /*clCreateKernelsInProgram(clDict->program,NumberOfKernels,clDict->kernels,&numkernels);
    printf("numkernels: %d\n",numkernels);*/
    for(i = 0; i < NumberOfKernels; ++i) {
        printf("%s\n",KERNELNAMES[i]);
        err = initKernel(clDict, KERNELNAMES[i], i);
        if(err != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    for(i=0; i<NumberOfKernels; ++i) {
        char name[512];
        cl_uint numargs;
        clGetKernelInfo(clDict->kernels[i],CL_KERNEL_FUNCTION_NAME,sizeof(name),name,
                        NULL);

        printf("info: %s\n",name);
        clGetKernelInfo(clDict->kernels[i],CL_KERNEL_NUM_ARGS,sizeof(numargs),&numargs,
                        NULL);
        printf("info: %d\n",numargs);
    }

    return EXIT_SUCCESS;
}

void createCLBuffer(CLDict *clDict, enum BUFFER buffer, size_t size, cl_mem_flags type){
    cl_int err;
    clDict->buffers[buffer] = clCreateBuffer(clDict->context,  type,
                                          size,NULL, &err);
    handleCLErr(err, clDict,"Error: Failed create buffer!");
}

void createCLBuffers(CLDict *clDict)
{
    createCLBuffer(clDict,QCL,sizeof(double)*QSIZE,CL_MEM_READ_WRITE);
    createCLBuffer(clDict,PCL,sizeof(double)*PSIZE,CL_MEM_READ_WRITE);
    createCLBuffer(clDict,LOGPCL,sizeof(double)*PSIZE,CL_MEM_READ_WRITE);

    if (FREQSCORR) {
        createCLBuffer(clDict,FSTCL,sizeof(double)*MAXPOPS,CL_MEM_READ_WRITE);
        createCLBuffer(clDict,EPSILONCL,sizeof(double)*NUMLOCI*MAXALLELES,CL_MEM_READ_WRITE);
    } else {
        createCLBuffer(clDict,LAMBDACL,sizeof(double)*MAXPOPS,CL_MEM_READ_WRITE);
    }

    createCLBuffer(clDict,ZCL,sizeof(int)*ZSIZE,CL_MEM_READ_WRITE);
    createCLBuffer(clDict,GENOCL,sizeof(int)*GENOSIZE,CL_MEM_READ_WRITE);


    if (RECESSIVEALLELES) {
        createCLBuffer(clDict,PREGENOCL,sizeof(int)*GENOSIZE,CL_MEM_READ_WRITE);
        createCLBuffer(clDict,RECESSIVECL,sizeof(int)*NUMLOCI,CL_MEM_READ_WRITE);

    }
    createCLBuffer(clDict,NUMALLELESCL,sizeof(int)*NUMLOCI,CL_MEM_READ_WRITE);

    if (PFROMPOPFLAGONLY || USEPOPINFO) {
        createCLBuffer(clDict,POPFLAGCL,sizeof(int)*NUMINDS,CL_MEM_READ_WRITE);
    }

    createCLBuffer(clDict,NUMAFROMPOPSCL,sizeof(int)*NUMLOCI*MAXPOPS*MAXALLELES,CL_MEM_READ_WRITE);

    createCLBuffer(clDict,NUMLOCIPOPSCL,sizeof(int)*NUMINDS*MAXPOPS,CL_MEM_READ_WRITE);

    createCLBuffer(clDict,LOGDIFFSCL,sizeof(double)*NUMINDS,CL_MEM_READ_WRITE);

    /*createCLBuffer(clDict,LOGTERMSCL,sizeof(double)*NUMINDS*NUMLOCI,CL_MEM_READ_WRITE);*/
    createCLBuffer(clDict,ALPHACL,sizeof(double)*MAXPOPS,CL_MEM_READ_WRITE);
    createCLBuffer(clDict,TESTQCL,sizeof(double)*QSIZE,CL_MEM_READ_WRITE);
    createCLBuffer(clDict,RANDCL,sizeof(double)*RANDSIZE,CL_MEM_READ_WRITE);
    createCLBuffer(clDict,ERRORCL,sizeof(int)*2,CL_MEM_READ_WRITE);

}


/*
 * Inits the dict, but the program must still be compiled.
 */
int InitCLDict(CLDict *clDictToInit)
{
    cl_kernel *kernels;
    cl_mem *buffers;
    size_t *locals;
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
    cl_context context;
    cl_device_id device_id;
    cl_command_queue commands;
    int DEVICETYPE;
    int err;
    int compileret;
    char options[1024];
    DEVICETYPE =  USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    err = clGetDeviceIDs(platform_id, DEVICETYPE, 1, &device_id, &ret_num_devices);
    if (err != CL_SUCCESS) {
        printf("retval %d\n",(int) ret);
        switch(err) {
        case CL_INVALID_PLATFORM:
            printf("invalid platform!");
            break;
        case CL_INVALID_VALUE:
            printf("invalid value");
            break;
        case CL_DEVICE_NOT_FOUND:
            printf("device not found");
            break;
        case CL_INVALID_DEVICE_TYPE:
            if(USEGPU) {
                printf("invalid device: GPU\n");
            } else {
                printf("invalid device: CPU\n");
            }
            break;
        }
        printf("Error: Failed to create a device group!\n");
        return EXIT_FAILURE;
    }

    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    if (!context) {
        printf("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }

    commands = clCreateCommandQueue(context, device_id, 0, &err);
    if (!commands) {
        printf("Error: Failed to create a command commands!\n");
        return EXIT_FAILURE;
    }



    kernels = calloc(NumberOfKernels, sizeof(cl_kernel));
    buffers = calloc(NumberOfBuffers, sizeof(cl_mem));
    locals = calloc(NumberOfKernels, sizeof(size_t));
    clDictToInit->kernels = kernels;
    clDictToInit->buffers = buffers;
    clDictToInit->locals = locals;
    clDictToInit->platform_id = platform_id;
    clDictToInit->ret_num_devices = ret_num_devices;
    clDictToInit->ret_num_platforms = ret_num_platforms;
    clDictToInit->device_id = device_id;
    clDictToInit->context = context;
    clDictToInit->commands = commands;

    /* compile OpenCL kernels */
    /*Define the constants in the kernels */
    /*
     * These are known and global when we compile the kernels
     * so we pass them over to the kernel compiler
     * to be inserted into the kernel code
     */



    sprintf(options,"-Werror -D UNASSIGNED=%d  -D MAXPOPS=%d -D MISSING=%d \
                     -D MAXALLELES=%d -D NUMLOCI=%d  -D LINES=%d    \
                     -D NUMINDS=%d -D MAXRANDOM=%d  -D USEPOPINFO=%d    \
                     -D LOCPRIOR=%d  -D NOTAMBIGUOUS=%d  -D NUMLOCATIONS=%d    \
                     -D PFROMPOPFLAGONLY=%d -D FREQSCORR=%d -D blockSize=64\
                     -D DEBUGCOMPARE=%d "
            , UNASSIGNED, MAXPOPS, MISSING
            , MAXALLELES, NUMLOCI, LINES
            , NUMINDS, MAXRANDOM, USEPOPINFO
            , LOCPRIOR, NOTAMBIGUOUS, NUMLOCATIONS
            , PFROMPOPFLAGONLY,FREQSCORR,DEBUGCOMPARE);

    printf("%s\n",options);
    compileret = CompileKernels(clDictToInit,options);

    if(compileret != EXIT_SUCCESS) {
        printf("Kernels failed to compile!\n");
        return EXIT_FAILURE;
    } else {
        printf("Kernels compiled!\n");
    }
    createCLBuffers(clDictToInit);
    printf("Buffers created!\n");
    setKernelArgs(clDictToInit);
    printf("Kernel args set!\n");
    getLocals(clDictToInit);
    printf("Work group info fetched!\n");
    return EXIT_SUCCESS;
}

int dimLoc(int * dims, int * dimMaxs, int numDims)
{
    int loc = 0;
    int i, j;
    for(i = 0; i < numDims; ++i) {
        int dimProd =1;
        for(j = i+1; j < numDims; j++) {
            dimProd *= dimMaxs[j];
        }
        loc += dims[i]*dimProd;
    }
    return loc;
}
/*
 * Copies the entire last dimension over into localarr
 */
void copyToLocal( double * globalArr, double *localArr,
                  int * dims, int * dimMaxs, int numDims)
{
    int i;
    int origLastDim = dims[numDims-1];
    for(i = 0; i < dimMaxs[numDims-1]; ++i) {
        dims[numDims-1] = i;
        localArr[i] = globalArr[dimLoc(dims,dimMaxs,numDims)];
    }
    dims[numDims-1] = origLastDim;
}

/*
 * Reads the buffer fource from the gpu to the array dest
 */
void readBuffer(CLDict *clDict, void * dest, size_t size, enum BUFFER source,
                char *name)
{
    cl_int err;
    char msg[120];
    err = clEnqueueReadBuffer(clDict->commands, clDict->buffers[source], CL_TRUE,
                              0,
                              size, dest, 0, NULL, NULL );
    strcpy(msg,"Failed to read buffer: ");
    strcat(msg,name);
    strcat(msg,"!\n");
    handleCLErr(err, clDict,msg);
    err = clFinish(clDict->commands);
    handleCLErr(err, clDict,"clFinish error!\n");
}

/*
 * Writes the array source to the buffer dest on the GPU
 */

void writeBuffer(CLDict *clDict, void * source, size_t size,
                 enum BUFFER dest, char *name)
{
    cl_int err;
    char msg[120];
    err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[dest], CL_TRUE,
                               0,
                               size, source, 0, NULL, NULL );
    strcpy(msg,"Failed to write buffer: ");
    strcat(msg,name);
    strcat(msg,"!\n");
    handleCLErr(err, clDict,msg);
    err = clFinish(clDict->commands);
    handleCLErr(err, clDict,"clFinish error!\n");
}

void runKernel(CLDict *clDict, enum KERNEL kernel, int numdims, size_t *dims,
               char *name)
{
    cl_int err;
    char msg[120];
    err = clEnqueueNDRangeKernel(clDict->commands, clDict->kernels[kernel],
                                 numdims, NULL, dims, NULL, 0, NULL, NULL);
    strcpy(msg,"Failed to run kernel: ");
    strcat(msg,name);
    strcat(msg,"!\n");
    handleCLErr(err, clDict,msg);
    err = clFinish(clDict->commands);
    handleCLErr(err, clDict,"clFinish error!\n");
}



/*int main(){
    CLDict *clDict = NULL;
    int ret;
    char *options = "-Werror -D UNASSIGNED=-1  -D MAXPOPS=2 -D MISSING=-999 \
                     -D MAXALLELES=15 -D NUMLOCI=15  -D LINES=2    \
                     -D NUMINDS=20 -D MAXRANDOM=2  -D USEPOPINFO=1    \
                     -D LOCPRIOR=1  -D NOTAMBIGUOUS=-1  -D NUMLOCATIONS=5    \
                     -D PFROMPOPFLAGONLY=0 -D blockSize=64";

    clDict = malloc(sizeof *clDict);

    InitCLDict(clDict);
    ret = CompileKernels(clDict,options);
    printf("return code: %d\n",ret);
    ReleaseCLDict(clDict);
    return ret;
}*/
