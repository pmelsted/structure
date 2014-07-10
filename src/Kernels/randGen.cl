#ifndef RANDGEN
#define RANDGEN
#define MAXRANDVAL 4294967296

#define RANDSPERSTREAM 2147483648
#define PI 3.14159265359

mwc64x_state_t getRandGen(__global uint *randGens, int id){
    /*mwc64x_state_t rng = { randGens[id*2] };*/

    /*mwc64x_state_t rng;
    int offset = randGen*2;
    rng.x = randGens[id*2];
    rng.c = randGens[id*2+1];*/

    /* C trickery */
    return  (mwc64x_state_t) { randGens[id*2] };
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


__kernel void InitRandGens( __global uint *randGens, const int baseOffset)
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

void initRndDiscState(RndDiscState *state, __global uint * randGens, int id)
{
    state->randGens = randGens;
    state->randGenId = id;
    state->rng = getRandGen(randGens,id);
}

double rndDisc(RndDiscState * state)
{
    double val = uintToUnit(MWC64X_NextUint(&(state->rng)));
    return val;
}

uint rndUInt(RndDiscState * state)
{
    uint val = MWC64X_NextUint(&(state->rng));
    return val;
}

double2 BoxMuller(RndDiscState *state)
{
  double u0=rndDisc(state), u1=rndDisc(state);
  double r=sqrt(-2*log(u0));
  double theta=2*PI*u1;
  return (double2) (r*sin(theta),r*cos(theta));
}

float2 BoxMullerF(RndDiscState *state)
{
  float u0=rndDisc(state), u1=rndDisc(state);
  float r=sqrt(-2*log(u0));
  float theta=2*PI*u1;
  return (float2) (r*sin(theta),r*cos(theta));
}



void saveRndDiscState(RndDiscState *state){
    saveRandGen(state->randGens,state->randGenId,state->rng);
}
/*#endif*/

#endif

