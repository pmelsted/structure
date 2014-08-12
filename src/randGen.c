#include "randGen.h"
#include "debug.h"
#include <stdio.h>
#include <stdlib.h>
#include "Kernels.h"
#include "ran.h"

void initRndDiscState(RndDiscState *state, float * randomArr, int maxrandom)
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


float rndDisc(RndDiscState * state)
{
    float val;

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

void FillArrayWithRandomCL(CLDict *clDict,float *randomArr, int numrands){
    size_t global[1];
    global[0] = numrands;
    setKernelArgExplicit(clDict,FillArrayWRandomKernel,sizeof(int),&numrands,2);
    runKernel(clDict,FillArrayWRandomKernel,1,global,"FillArrayWRandom");
    /*readBuffer(clDict,randomArr, sizeof(float) * numrands,RANDCL,*/
                /*"randomArr");*/
}

void FillArrayWithRandom(float *random, int n)
{
    /*Fills an array with floating random numbers in [0,1)*/
    int i;
    float d;
    for(i = 0; i < n; i++) {
        d = (float) rand() / ((float) RAND_MAX + 1);
        random[i] = d;
    }
}
