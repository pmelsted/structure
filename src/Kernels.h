#include "KernelDefs.h"
extern int InitCLDict(CLDict *clDictToInit);
extern void ReleaseCLDict(CLDict *clDict);
extern int CompileKernels(CLDict *clDict, char * names[], char *vals[],int numVals);
