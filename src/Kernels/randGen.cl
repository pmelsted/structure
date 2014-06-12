typedef struct RndDiscState {
    __global double *randomArr;
    int randomValsTaken;
    int maxrandom;
    int baseOffset;
} RndDiscState;


void initRndDiscState(RndDiscState *state,__global double * randomArr,
                     int maxrandom)
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
    return val;
}
