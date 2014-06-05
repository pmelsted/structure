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
#define RandPos(ind,loc) ind*NUMLOCI+loc

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
              Z[ZPos (ind, line, loc)] = PickAnOptionDiscrete (MAXPOPS, sum, Cutoffs,randArr[RandPos(ind,loc)]);
           }
       }
   }
}


#if LINES == 2
__kernel void UpdateGeno (
   __global double* Q, /* input */
   __global double* P,  /* input */
   __global int* PreGeno, /* input */
   __global int* Recessive, /* input */
   __global int* NumAlleles, /* input */
   __global int* Z, /* input */
   __global double* randArr, /*random numbers*/
   __global int* Geno/* output */
   )
{
   int dom;
   double AlleleProbs[4];
   double Sum;

   int ind = get_global_id(0);
   int loc = get_global_id(1); /* is this correct? */

   if(ind < NUMINDS && loc < NUMLOCI){
       if (PreGeno[GenPos (ind, 0, loc)] != MISSING
           && PreGeno[GenPos (ind, 1, loc)] != MISSING) {
            for (dom = 0; dom < 4; dom++) {
              AlleleProbs[dom] = 0.0;
            }

            #if RECESSIVEALLELES
            /* known at build time */
            if (Recessive[loc] != MISSING    /* bug fixed 05072007 */
                && PreGeno[GenPos (ind, 0, loc)] != Recessive[loc]) {
              AlleleProbs[0] =
                  P[PPos(loc, Z[ZPos(ind, 0, loc)], Recessive[loc])] *
                  P[PPos(loc, Z[ZPos(ind, 1, loc)], PreGeno[GenPos(ind, 0, loc)])];
              AlleleProbs[1] =
                  P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                  P[PPos(loc, Z[ZPos (ind, 1, loc)], Recessive[loc])];
            }
            #endif

            AlleleProbs[2] =
                P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                P[PPos(loc, Z[ZPos (ind, 1, loc)], PreGeno[GenPos (ind, 1, loc)])];

            Sum = AlleleProbs[0] + AlleleProbs[1] + AlleleProbs[2];
            dom = PickAnOptionDiscrete (3, Sum, AlleleProbs,randArr[RandPos(ind,loc)]);

            if (dom == 0) {
              Geno[GenPos (ind, 0, loc)] = Recessive[loc];
            } else {
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            }
            if (dom == 1){
                Geno[GenPos (ind, 1, loc)] = Recessive[loc];
            } else {
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
            }

      }
   }
}
#else
__kernel void UpdateGeno (
   __global double* Q, /* input */
   __global double* P,  /* input */
   __global int* PreGeno, /* input */
   __global int* Recessive, /* input */
   __global int* NumAlleles, /* input */
   __global int* Z, /* input */
   __global double* randArr, /*random numbers*/
   __global int* Geno/* output */
   )
{
    /*undefined*/
}
#endif
