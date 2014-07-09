#include "randGen.h"
#include "debug.h"
#include <stdio.h>
#include <stdlib.h>
#include "Kernels.h"
#include "ran.h"

void initRndDiscState(RndDiscState *state, double * randomArr, int maxrandom)
{
    state->randomArr = randomArr;
    state->maxrandom = maxrandom;
    state->baseOffset = 0;
    state->randomValsTaken = 0;
}

void rndDiscStateReset(RndDiscState *state, int baseOffset)
{
    state->baseOffset = baseOffset;
    state->randomValsTaken = 0;
}


double rndDisc(RndDiscState * state)
{
    double val;

    val = state->randomArr[state->baseOffset + state->randomValsTaken];
    state->randomValsTaken++;
    if (state->randomValsTaken > state->maxrandom) {
        printf("Too many random vals taken!\n");
        printf("Random values taken %d\n", state->randomValsTaken -1);
        print_trace();
        exit(1);
    }
    return val;
}

void FillArrayWithRandomCL(CLDict *clDict,double *randomArr, int numrands){
    size_t global[1];
    global[0] = numrands/16;
    setKernelArgExplicit(clDict,FillArrayWRandomKernel,sizeof(int),&numrands,2);
    runKernel(clDict,FillArrayWRandomKernel,1,global,"FillArrayWRandom");
    /*readBuffer(clDict,randomArr, sizeof(double) * numrands,RANDCL,*/
                /*"randomArr");*/
}

void FillArrayWithRandom(double *random, int n)
{
    /*Fills an array with floating random numbers in [0,1)*/
    int i;
    double d;
    for(i = 0; i < n; i++) {
        d = (double) rand() / ((double) RAND_MAX + 1);
        random[i] = d;
    }
}
