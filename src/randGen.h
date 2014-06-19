#ifndef randGen
#define randGen
typedef struct RndDiscState {
    double *randomArr;
    int randomValsTaken;
    int maxrandom;
    int baseOffset;
} RndDiscState;

extern void initRndDiscState(RndDiscState *state,
                             double * randomArr, int maxrandom);

extern double rndDisc(RndDiscState * state);
extern void rndDiscStateReset(RndDiscState *state,
                              int baseOffset);
#endif

