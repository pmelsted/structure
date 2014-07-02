#include "Kernels/rand.cl"
__kernel void UpdateZ (
    __global double* Q, /* input */
    __global double* P,  /* input */
    __global int* Geno,/* input */
    __global uint* randGens, /*random numbers*/
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


    if(ind < NUMINDS && loc < NUMLOCI) {
        initRndDiscState(randState,randGens,ind*NUMLOCI+loc);
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
                Z[ZPos (ind, line, loc)] = PickAnOptionDiscrete (MAXPOPS, sum, Cutoffs,
                                           randState);
            }
        }
        saveRndDiscState(randState);
    }
}

