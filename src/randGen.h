#ifndef randGen
#define randGen
#include "Kernels.h"
typedef struct RndDiscState {
    double *randomArr;
    int randomValsTaken;
    int maxrandom;
    int baseOffset;
} RndDiscState;

extern void initRndDiscState(RndDiscState *state, double * randomArr,
                             int maxrandom);

extern double rndDisc(RndDiscState * state);
extern void rndDiscStateReset(RndDiscState *state, int baseOffset);
extern void FillArrayWithRandomCL(CLDict *clDict,double *randomArr, int numrands);
extern void FillArrayWithRandom(double *random, int n);
#endif

