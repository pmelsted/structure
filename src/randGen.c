#include "randGen.h"
#include <stdio.h>
#include <stdlib.h>

void initRndDiscState(RndDiscState *state, double * randomArr,
                     int maxrandom, int baseOffset)
{
    state->randomArr = randomArr;
    state->maxrandom = maxrandom;
    state->baseOffset = baseOffset;
    state->randomValsTaken = 0;
}

double rndDisc(RndDiscState * state)
{
    double val;
    val = state->randomArr[state->baseOffset + state->randomValsTaken];
    state->randomValsTaken++;
    if (state->randomValsTaken > state->maxrandom){
        printf("Too many random vals taken!\n");
        printf("Random values taken %d\n", state->randomValsTaken -1);
        exit(1);
    }
    return val;
}
