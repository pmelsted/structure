#include "randGen.h"
#include "debug.h"
#include <stdio.h>
#include <stdlib.h>

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
