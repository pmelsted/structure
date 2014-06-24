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
#endif

