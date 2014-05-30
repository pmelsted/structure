#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>
#include "KernelDefs.h"

#define MAX_SOURCE_SIZE (0x100000)
#define USEGPU 1



/*
 * Inits the dict, but the program must still be compiled.
 */
int InitCLDict(CLDict *clDictToInit){
    cl_kernel kernels[NumberOfKernels];    
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
    cl_context context;
    cl_device_id device_id; 
    cl_command_queue commands;
    int DEVICETYPE;
    int err;
    DEVICETYPE =  USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
    err = clGetDeviceIDs(platform_id, DEVICETYPE, 1, &device_id, &ret_num_devices);
    if (err != CL_SUCCESS)
    {
        printf("retval %d\n",(int) ret);
    switch(err){
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
            if(USEGPU){
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
    if (!context)
    {
        printf("Error: Failed to create a compute context!\n");
        return EXIT_FAILURE;
    }

    commands = clCreateCommandQueue(context, device_id, 0, &err);
    if (!commands)
    {
        printf("Error: Failed to create a command commands!\n");
        return EXIT_FAILURE;
    }




    clDictToInit->numKernelsInDict = 0;
    clDictToInit->kernels = kernels;
    clDictToInit->platform_id = platform_id;
    clDictToInit->ret_num_devices = ret_num_devices;
    clDictToInit->ret_num_platforms = ret_num_platforms;
    clDictToInit->device_id = device_id;
    clDictToInit->context = context;
    clDictToInit->commands = commands;
    return EXIT_SUCCESS;
}

void ReleaseCLDict(CLDict *clDict){
    int i;
    clReleaseProgram(clDict->program);
    for(i = 0; i < clDict->numKernelsInDict; ++i){
        clReleaseKernel(clDict->kernels[i]);
    } 
    clReleaseCommandQueue(clDict->commands);
    clReleaseContext(clDict->context);
    free(clDict);
}



char * searchReplace(char * string,  char *toReplace[], char *replacements[], int numReplacements){
    int i = 0;
    char *locOfToRep;
    char *toRep;
    char *rep;
    int lenToRep,lenStr,lenAfterLocRep;
    static char buffer[MAX_SOURCE_SIZE];
    strncpy(buffer,string,MAX_SOURCE_SIZE);
    for(i = 0; i < numReplacements; ++i){
        toRep = toReplace[i];
        rep = replacements[i];
        /*if str not in the string, skip it */
        if (!(locOfToRep = strstr(buffer,toRep))){
           printf("Warning: %s not found in string!\n",toRep);
           continue;
        }
        lenToRep = strlen(toRep); 
        lenStr = strlen(buffer); 
        lenAfterLocRep = strlen(locOfToRep); 
        /*Print the string upto the pointer, then the val, and then the rest of the string.*/
        sprintf(buffer, "%.*s%s%s", lenStr-lenAfterLocRep, buffer,rep,locOfToRep+lenToRep);
        /*Will break if buffer is longer than string, so we restrict ourselves to the longest length.*/
    }
    return buffer;
}

void preProcessSource(char * kernelSource, size_t *source_size, char *names[], char *vals[], int numVals){
   char * processedSource;
   processedSource = searchReplace(kernelSource, names, vals, numVals);
   strncpy(kernelSource,processedSource,MAX_SOURCE_SIZE);
   *source_size = strlen(kernelSource);
}



int initKernel(CLDict *clDict,char * kernelName, enum KERNEL kernelEnumVal){
    cl_kernel kernel;
    int err;
    kernel = clCreateKernel(clDict->program, kernelName, &err);
    if (!kernel || err != CL_SUCCESS)
    {
    	switch(err){
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
int CompileKernels(CLDict *clDict, char *names[],char *vals[], int numVals){
    FILE *fp;
    char *KernelSource;
    size_t source_size;
    
    int err;
    int i;
    cl_program program;

    char *kernelNames[NumberOfKernels] = {"UpdateZ"};
    /*cl_int ret;*/



    /* Load the source code containing the kernels*/
    fp = fopen("Kernels/Kernels.cl", "r");
    if (!fp) {
    	fprintf(stderr, "Failed to load kernel file Kernels.cl\n");
        return EXIT_FAILURE;
    }


    KernelSource = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(KernelSource, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    preProcessSource(KernelSource, &source_size, names,vals,numVals); 
    program = clCreateProgramWithSource(clDict->context, 1, (const char **) & KernelSource, NULL, &err);
    if (!program)
    {
        printf("Error: Failed to create compute program!\n");
        return EXIT_FAILURE;
    }


    free(KernelSource);


    err = clBuildProgram(program, 1, &clDict->device_id, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];

        printf("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, clDict->device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        printf("%d\n", (int) len);
        return EXIT_FAILURE;
    }

    clDict->program = program;
    for(i = 0; i < NumberOfKernels; ++i){
        err = initKernel(clDict, kernelNames[i], i);
        if(err != EXIT_SUCCESS){
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

/*
int main(int argc, char *argv[]){
    CLDict *clDict = NULL;
    int numVals;
    int ret;
    char *names[5] = {"%maxpops%", "%missing%", "%maxalleles%","%numloci%","%lines%"};
    char *vals[5] = {"2", "-999", "15","15","2"};
    clDict = malloc(sizeof *clDict);
    
    numVals = 5;
    InitCLDict(clDict);
    ret = CompileKernels(clDict,names,vals,numVals);
    printf("return code: %d\n",ret);
    ReleaseCLDict(clDict);
    return ret;
}*/
