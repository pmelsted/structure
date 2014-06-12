#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define UNASSIGNED -9

//These are inserted during PreProcessingDuringCompilation.
#define MISSING %missing%
#define MAXPOPS %maxpops%
#define MAXALLELES %maxalleles%
#define NUMLOCI %numloci%
#define LINES %lines%
#define NUMINDS %numinds%
#define MAXRANDOM %maxrandom%
#define NOTAMBIGUOUS %notambiguous%
#define USEPOPINFO %usepopinfo%
#define LOCPRIOR %locprior%
#define NUMLOCATIONS %numlocations%

#define GenPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))

#include "Kernels/util.cl"
#include "Kernels/KernelErrors.h"
#include "Kernels/UpdateZ.cl"
#include "Kernels/UpdateQ.cl"

