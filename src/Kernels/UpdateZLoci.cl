
#include "Kernels/util.cl"
#define UNASSIGNED -9

//These are inserted during PreProcessingDuringCompilation.
#define MISSING %missing%
#define MAXPOPS %maxpops%
#define MAXALLELES %maxalleles%
#define NUMLOCI %numloci%
#define LINES %lines%

#define GenPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))


__kernel void UpdateZLoci (                                                       
   __global float* Q, /* input */
   __global float* P,  /* input */                                           
   __global float* Geno,/* input */                                   
   __global float* Z, /* output */                                                   
   const unsigned int ind,
   const unsigned int line
   )                                      
{                                                                      
   int allele;
   int pop;
   float Cutoffs[MAXPOPS];
   float sum;
   int samplesPerStream = 100000;
   int baseOffset = 0; 
   mwc64x_state_t rng;
   //ulong samplesPerStream=n/get_global_size(0);
   MWC64X_SeedStreams(&rng, baseOffset, 2*samplesPerStream);
   
   int loc = get_global_id(0); /* is this correct? */

   allele = Geno[GenPos (ind,line,loc)];
   if(allele == MISSING){
       Z[ZPos(ind,line,loc)] = UNASSIGNED;
   } else {
       sum = 0.0;
       for (pop = 0; pop < MAXPOPS; pop++) {
         Cutoffs[pop] = Q[QPos (ind, pop)] * P[PPos (loc, pop, allele)];
         sum += Cutoffs[pop];
       }
       Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs,&rng);
   }
}                                                                     
