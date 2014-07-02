#ifndef RANDGEN
#define RANDGEN
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


mwc64x_state_t getRandGen(__global uint *randGens, int randGen){
    mwc64x_state_t rng;
    int offset = randGen*2;
    uint x = randGens[offset];
    uint c = randGens[offset+1];
    rng.x = x;
    rng.c = c;
    return rng;
}
void saveRandGen(__global uint *randGens, int randGen,mwc64x_state_t rng){
    uint x = rng.x;
    uint c = rng.c;
    int offset = randGen*2;
    randGens[offset]=x;
    randGens[offset+1]=c;

}
__kernel void InitRandGens(
        __global uint *randGens,
        const int baseOffset
        )
{
    int pos = get_global_id(0);
    mwc64x_state_t rng;
    MWC64X_SeedStreams(&rng,baseOffset,2^20);
    saveRandGen(randGens,pos,rng);
}
#endif
