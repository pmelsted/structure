#ifndef Kernels
#define Kernels
#include "KernelDefs.h"
extern int InitCLDict(CLDict *clDictToInit);
extern void ReleaseCLDict(CLDict *clDict);
extern int CompileKernels(CLDict *clDict, char * names[], char *vals[],
                          int numVals);
extern void printCLErr(cl_int err);
extern void PrintKernelError(int error);
extern void handleCLErr(cl_int err,CLDict *clDict,char * message);
extern void createCLBuffers(CLDict *clDict);
extern void copyToLocal( double * globalArr, double *localArr,
                         int * dims, int * dimMaxs, int numDims);
extern int dimLoc(int * dims, int * dimMaxs, int numDims);
extern void readBuffer(CLDict *clDict, void * dest, size_t size,
                       enum BUFFER buffer, char *name);
extern void readBuffers(CLDict *clDict, void * dest[], size_t size[], enum BUFFER source[],
                char *name[], int numbufferstoread);
extern void writeBuffer(CLDict *clDict, void * source, size_t size,
                        enum BUFFER buffer, char *name);
extern void runKernel(CLDict *clDict, enum KERNEL kernel, int numdims,
                      size_t *dims, char *name);
extern void runTask(CLDict *clDict, enum KERNEL kernel, char *name);

extern void setKernelArgExplicit(CLDict *clDict, enum KERNEL kernel,size_t size, void *arg, int argnum);
extern void setKernelArg(CLDict *clDict, enum KERNEL kernel, enum BUFFER buffer,int argnum);

extern void finishCommands(CLDict *clDict, char * name);

#endif


