#define MAXRANDVAL 4294967296

__kernel void Dirichlet(
        __global double *Parameters,
        __global double *randomArr,
        __global double *TestQ)
{
    int ind = get_global_id(0);
    RndDiscState randState[1];

    initRndDiscState(randState,randomArr,MAXRANDOM);
    rndDiscStateReset(randState,ind*MAXRANDOM);
    double GammaSample[MAXPOPS];

    int i = 0;
    double sum = 0.0;
    for(i = 0; i < MAXPOPS; i++){
        GammaSample[i] = RGammaDisc(Parameters[i],1,randState);
        sum += GammaSample[i];
    }
    for(i = 0; i < MAXPOPS; i++){
        TestQ[i+ind*MAXPOPS] = GammaSample[i]/sum;
    }
}

double uintToUnit(uint rndint){
   return (double) rndint / MAXRANDVAL;
}

__kernel void FillArrayWRandom(
        __global double *randomArr,
        const int length,
        const int baseOffset
        )
{
    int pos = get_global_id(0);
    uint i;
    double val;
    ulong samplesPerstream = length/get_global_size(0);
    uint offset = pos*samplesPerstream;
    if (offset < length){
        mwc64x_state_t rng;
        MWC64X_SeedStreams(&rng,baseOffset,samplesPerstream);
        for(i = 0; i < samplesPerstream && i+offset < length; i++){
            randomArr[offset+i] = uintToUnit(MWC64X_NextUint(&rng));
        }
    }
}

