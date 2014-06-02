#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#include "Kernels/util.cl"
#define UNASSIGNED -9

//These are inserted during PreProcessingDuringCompilation.
#define MISSING %missing%
#define MAXPOPS %maxpops%
#define MAXALLELES %maxalleles%
#define NUMLOCI %numloci%
#define LINES %lines%
#define NUMINDS %numinds%

#define GenPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))

__kernel void UpdateZ (                                                       
   __global double* Q, /* input */
   __global double* P,  /* input */                                           
   __global int* Geno,/* input */
   __global double* randArr, /*random numbers*/                                   
   __global int* Z /* output */                                                   
   )                                      
{                                                                      
   int allele;
   int pop;
   int line;
   double Cutoffs[MAXPOPS];
   double sum;
   
   int ind = get_global_id(0);
   int loc = get_global_id(1); /* is this correct? */

   if(ind < NUMINDS && loc < NUMLOCI){
       for (line = 0; line < LINES; line++) {
           allele = Geno[GenPos (ind,line,loc)];
           if (allele == MISSING) {   /*Missing Data */
             Z[ZPos (ind, line, loc)] = UNASSIGNED;
           } else {
             /*Data present */
             sum = 0.0;    /*compute prob of each allele being from each pop */
              for (pop = 0; pop < MAXPOPS; pop++) {
                Cutoffs[pop] = Q[QPos (ind, pop)] * P[PPos (loc, pop, allele)];
                sum += Cutoffs[pop];
              }
              Z[ZPos (ind, line, loc)] = PickAnOptionDiscrete (MAXPOPS, sum, Cutoffs,randArr[ind*NUMLOCI+loc]);
           }
       }
   }
}                                                                    
