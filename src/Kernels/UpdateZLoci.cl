#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))

#include "Kernels/util.cl"
#define UNASSIGNED -9

//These are inserted during PreProcessingDuringCompilation.
#define MISSING %missing%
#define MAXPOPS %maxpops%

__kernel void UpdateZLoci (                                                       
   __global float* Q, /* input */
   __global float* P,  /* input */                                           
   __global float* Geno,/* input */                                   
   __global float* Z, /* output */                                                   
   const unsigned int NUMLOCI,
   const unsigned int ind,
   const unsigned int line
   )                                      
{                                                                      
   int allele;
   float Cutoffs[MAXPOPS];
   int id = get_global_id(0);                                           
   /*
   allele = Geno[GenPos (ind,line,loc)]
   if(allelle == MISSING){
       Z[ZPos(ind,line,loc)] = UNASSIGNED;
   } else {
       sum = 0.0;
       for (pop = 0; pop < MAXPOPS; pop++) {
         Cutoffs[pop] = Q[QPos (ind, pop)] * P[PPos (loc, pop, allele)];
         sum += Cutoffs[pop];
       }
       Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs);
   }
  */
}                                                                     
