#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//These are inserted during PreProcessingDuringCompilation.

#define GenPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))
#define NumAFromPopPos(pop,allele) ((pop)*(MAXALLELES)+(allele))    /* NumAFromPop */

#include "Kernels/util.cl"
#include "Kernels/KernelErrors.h"
#include "Kernels/UpdateZ.cl"
#include "Kernels/UpdateP.cl"

