#include "Kernels/rand.cl"
__kernel void UpdateZ (
    __global float* Q, /* input */
    __global float* P,  /* input */
    __global int* Geno,/* input */
    __global uint* randGens, /*random numbers*/
    __global int* Z, /* output */
    __global int* error
)
{
    int allele;
    int pop;
    int line;
    float Cutoffs[MAXPOPS];
    float sum;

    int ind = get_global_id(0);
    RndDiscState randState[1];


    while(ind < NUMINDS){
        int loc = get_global_id(1);
        while (loc < NUMLOCI) {
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
            loc += get_global_size(1);
        }
        ind += get_global_size(0);
    }
}

