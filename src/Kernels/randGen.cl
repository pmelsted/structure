#ifndef RANDGEN
#define RANDGEN
#define MAXRANDVAL 4294967296

#define RANDSPERSTREAM 2147483648

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

double uintToUnit(uint rndint){
   return (double) rndint / MAXRANDVAL;
}


__kernel void InitRandGens(
        __global uint *randGens,
        const int baseOffset
        )
{
    int pos = get_global_id(0);
    mwc64x_state_t rng;
    MWC64X_SeedStreams(&rng,baseOffset,RANDSPERSTREAM);
    saveRandGen(randGens,pos,rng);
}

/*#if DEBUGCOMPARE
typedef struct RndDiscState {
    __global double *randomArr;
    int randomValsTaken;
    int baseOffset;
} RndDiscState;


void initRndDiscState(RndDiscState *state,__global double * randomArr, int offset)
{
    state->randomArr = randomArr;
    state->maxrandom = maxrandom;
    state->baseOffset = offset;
    state->randomValsTaken = 0;
}

double rndDisc(RndDiscState * state)
{
    double val;
    val = state->randomArr[state->baseOffset + state->randomValsTaken];
    state->randomValsTaken++;
    return val;
}
#else*/
typedef struct RndDiscState {
    __global uint *randGens;
    int randGenId;
    mwc64x_state_t rng;
} RndDiscState;

void initRndDiscState(RndDiscState *state,
                      __global uint * randGens,
                      int offset)
{
    state->randGens = randGens;
    state->randGenId = offset;
    state->rng = getRandGen(randGens,offset);
}

double rndDisc(RndDiscState * state)
{
    double val = uintToUnit(MWC64X_NextUint(&(state->rng)));
    return val;
}

void saveRndDiscState(RndDiscState *state){
    saveRandGen(state->randGens,state->randGenId,state->rng);
}
/*#endif*/

#endif

