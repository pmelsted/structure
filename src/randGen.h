#ifndef randGen
#define randGen
#include "Kernels.h"
typedef struct RndDiscState {
    float *randomArr;
    int randomValsTaken;
    int maxrandom;
    int baseOffset;
} RndDiscState;

extern void initRndDiscState(RndDiscState *state, float * randomArr,
                             int maxrandom);

extern float rndDisc(RndDiscState * state);
extern void rndDiscStateReset(RndDiscState *state, int baseOffset);
extern void FillArrayWithRandomCL(CLDict *clDict,float *randomArr, int numrands);
extern void FillArrayWithRandom(float *random, int n);
#endif

