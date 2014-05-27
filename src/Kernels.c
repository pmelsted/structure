#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <CL/cl.h>

#define MAXNUMBEROFKERNELS 1
#define MAX_SOURCE_SIZE (0x100000)
#define USEGPU 1


typedef struct KernelDict {
    int numKernelsInDict;
    char **KernelNames;
    cl_program *Kernels;    
} KernelDict;


/*
 * If we've already compiled this kernel,
 * return the index in CompiledKernels in which it is located.
 */
int alreadyCompiled(char * kernelFilename, KernelDict * alreadyCompiledKernels){
    int i, comp;
    for(i = 0; i < (*alreadyCompiledKernels).numKernelsInDict && i < MAXNUMBEROFKERNELS; i++){
      comp = strcmp(kernelFilename,(*alreadyCompiledKernels).KernelNames[i]); 
      if(comp == 0){
        return i;     
      }
    }
    return -1;
}


int CompileProgram(char * kernelFilename, KernelDict *alreadyCompiledKernels){
    int indIfAlreadyCompiled;
    indIfAlreadyCompiled = alreadyCompiled(kernelFilename, alreadyCompiledKernels);
    if(indIfAlreadyCompiled >= 0){
        return indIfAlreadyCompiled;        
    }

    FILE *fp;
    char *KernelSource; //Pointer to the kernel source
    size_t source_size; //Size of the kernel file (not used)
    
    int err;
    cl_program program;
    cl_platform_id platform_id;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;
    cl_int ret;
    cl_context context;
    cl_device_id device_id; 
    int DEVICETYPE;

    /* Load the source code containing the kernel*/
    fp = fopen(kernelFilename, "r");
    if (!fp) {
    	fprintf(stderr, "Failed to load kernel %s.\n", kernelFilename);
        exit(EXIT_FAILURE);
    }

    platform_id = NULL;
    KernelSource = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread(KernelSource, 1, MAX_SOURCE_SIZE, fp);
    fclose(fp);

    ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
	
 
    DEVICETYPE =  USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;
    err = clGetDeviceIDs(platform_id, DEVICETYPE, 1, &device_id, &ret_num_devices);
    if (err != CL_SUCCESS)
    {
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
        exit(EXIT_FAILURE);
    }
  
    // Create a compute context 
    //
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    if (!context)
    {
        printf("Error: Failed to create a compute context!\n");
        exit(EXIT_FAILURE);
    }

    // Create the compute program from the source buffer
    //
    program = clCreateProgramWithSource(context, 1, (const char **) & KernelSource, NULL, &err);
    if (!program)
    {
        printf("Error: Failed to create compute program!\n");
        exit(EXIT_FAILURE);
    }

    // Build the program executable
    //
    err = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];

        printf("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(EXIT_FAILURE);
    }

    int indOfKernel;

    indOfKernel = (*alreadyCompiledKernels).numKernelsInDict++;
    (*alreadyCompiledKernels).KernelNames[indOfKernel] = kernelFilename;
    (*alreadyCompiledKernels).Kernels[indOfKernel] = program;
    return indOfKernel;
}

void InitKernelDict(KernelDict *kernelDictToInit){
   char *KernelNames[MAXNUMBEROFKERNELS];
   cl_program Kernels[MAXNUMBEROFKERNELS];    

   (*kernelDictToInit).numKernelsInDict = 0;
   (*kernelDictToInit).KernelNames = KernelNames;
   (*kernelDictToInit).Kernels = Kernels;
}


int main(int argc, char *argv[]){
    int i;
    KernelDict alreadyCompiledKernels;
    InitKernelDict(&alreadyCompiledKernels);
    int a = CompileProgram("UpdateZLoci.cl", &alreadyCompiledKernels); 
    int b = CompileProgram("UpdateZLoci.cl", &alreadyCompiledKernels); 
    printf("Compiled!\n");
    printf("a %d, b %d\n",a,b);
}
