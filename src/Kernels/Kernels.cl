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
//#define NOTAMBIGUOUS %notambiguous%
#define NOTAMBIGUOUS -1

#define GenPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define ZPos(ind,line,loc) ((ind)*(LINES)*(NUMLOCI)+(line)*(NUMLOCI)+(loc))
#define PPos(loc,pop,allele) ((loc)*(MAXPOPS)*(MAXALLELES)+(pop)*(MAXALLELES)+(allele))
#define QPos(ind,pop) ((ind)*(MAXPOPS)+(pop))

#include "Kernels/util.cl"
#include "Kernels/KernelErrors.h"

__kernel void UpdateZ (
   __global double* Q, /* input */
   __global double* P,  /* input */
   __global int* Geno,/* input */
   __global double* randArr, /*random numbers*/
   __global int* Z, /* output */
   __global int* error
   )
{
   int allele;
   int pop;
   int line;
   double Cutoffs[MAXPOPS];
   double sum;

   int ind = get_global_id(0);
   int loc = get_global_id(1); /* is this correct? */
   /* we can't malloc in opencl, but we need only space for one */
   RndDiscState randState[1];

   initRndDiscState(randState,randArr,LINES);
   if(ind < NUMINDS && loc < NUMLOCI){
       rndDiscStateReset(randState,ind*NUMLOCI*LINES + loc*LINES);
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
              Z[ZPos (ind, line, loc)] = PickAnOptionDiscrete (MAXPOPS, sum, Cutoffs,randState);
           }
       }
       if (randState->randomValsTaken > randState->maxrandom){
           *error = KERNEL_OUT_OF_BOUNDS;
       }
   }
}


/*[>
 *  untested!
 <]
#if LINES == 2
__kernel void UpdateGeno (
   __global double* Q, [> input <]
   __global double* P,  [> input <]
   __global int* PreGeno, [> input <]
   __global int* Recessive, [> input <]
   __global int* NumAlleles, [> input <]
   __global int* Z, [> input <]
   __global double* randArr, [>random numbers<]
   __global int* Geno[> output <]
   )
{
    int dom;
    double AlleleProbs[4];
    double Sum;

    int ind = get_global_id(0);
    int loc = get_global_id(1);

    if(ind < NUMINDS && loc < NUMLOCI){
        if (PreGeno[GenPos (ind, 0, loc)] != MISSING
            && PreGeno[GenPos (ind, 1, loc)] != MISSING) {

            for (dom = 0; dom < 4; dom++) {
                AlleleProbs[dom] = 0.0;
            }
            [>
             * this will always be true, as UpdateGeno is not run if
             * there arent recessive alleles
             <]
            if (Recessive[loc] != MISSING    [> bug fixed 05072007 <]
                && PreGeno[GenPos (ind, 0, loc)] != Recessive[loc]) {
                AlleleProbs[0] =
                    P[PPos(loc, Z[ZPos(ind, 0, loc)], Recessive[loc])] *
                    P[PPos(loc, Z[ZPos(ind, 1, loc)], PreGeno[GenPos(ind, 0, loc)])];
                AlleleProbs[1] =
                    P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                    P[PPos(loc, Z[ZPos (ind, 1, loc)], Recessive[loc])];
            }
            AlleleProbs[2] =
                P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                P[PPos(loc, Z[ZPos (ind, 1, loc)], PreGeno[GenPos (ind, 1, loc)])];

            Sum = AlleleProbs[0] + AlleleProbs[1] + AlleleProbs[2];
            dom = PickAnOptionDiscrete (3, Sum, AlleleProbs,randArr[RandPos(ind,0,loc)]);

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
   __global double* Q, [> input <]
   __global double* P,  [> input <]
   __global int* PreGeno, [> input <]
   __global int* Recessive, [> input <]
   __global int* NumAlleles, [> input <]
   __global int* Z, [> input <]
   __global double* randArr, [>random numbers<]
           const int MAXREJECTIONS,
   __global int* Geno,[> output <]
   )
{

    [>
     * not yet fully implemented
     <]

    [>int ind = get_global_id(0);
    int loc = get_global_id(1);
    double AlleleProbs[MAXALLELES];
    int AlleleUsed[MAXALLELES];
    int AllelePresent[MAXALLELES];
    int alleleP;
    int allelecount;
    int notmissingcount;
    int toggle;
    int rejectioncount;

    if(ind < NUMINDS && loc < NUMLOCI){
        if (Recessive[loc]==NOTAMBIGUOUS) {
            for (line=0;line<LINES;line++) {
                Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
            }
        } else {
            allelecount=0;
            notmissingcount=0;

            for (allele = 0; allele < NumAlleles[loc]; allele++) {
                AllelePresent[allele] = 0;
            }

            for (line = 0; line < LINES; line++) {
                if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
                    notmissingcount+=1;
                    alleleP = AllelePresent[PreGeno[GenPos (ind, line, loc)]];
                    [> 0 if already present, otherwise 1 <]
                    AllelePresent[PreGeno[GenPos (ind, line, loc)]] += (1-allelP)*1;
                    allelecount+=(1 - alleleP);
              }
            }

            if (allelecount==notmissingcount) {  [> if number of alleles equal to number of slots then nothing to do <]
                for (line=0;line<LINES;line++) {
                    Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
                }
            } else {

            }
        }
    }<]
}
#endif*/
