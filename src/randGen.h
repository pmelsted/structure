#ifndef randGen
#define randGen
typedef struct RndDiscState {
    double *randomArr;
    int randomValsTaken;
    int maxrandom;
    int baseOffset;
} RndDiscState;

extern void initRndDiscState(RndDiscState *state, double * randomArr,
                     int maxrandom, int baseOffset);

extern double rndDisc(RndDiscState * state);
#endif
