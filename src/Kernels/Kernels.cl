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


__kernel void UpdateZ (                                                       
   __global float* Q, /* input */
   __global float* P,  /* input */                                           
   __global float* Geno,/* input */
   __global float* randNums, /*random numbers*/                                   
   __global float* Z, /* output */                                                   
   const unsigned int baseOffset
   )                                      
{                                                                      
   int allele;
   int pop;
   int line;
   float Cutoffs[MAXPOPS];
   float sum;
   mwc64x_state_t rng;
   //ulong samplesPerStream=n/get_global_size(0);
   ulong samplesPerStream= 10000;
   MWC64X_SeedStreams(&rng, baseOffset, 2*samplesPerStream);
   
   int ind = get_global_id(0);
   int loc = get_global_id(1); /* is this correct? */

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
          Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs,&rng);
       }
   }
}                                                                     
