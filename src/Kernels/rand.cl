#include "Kernels/randGen.cl"

__kernel void Dirichlet(
        __global float *Parameters,
        __global uint *randGens,
        __global float *TestQ)
{
    int ind = get_global_id(0);
    RndDiscState randState[1];

    while (ind < NUMINDS){
        initRndDiscState(randState,randGens,ind);
        float GammaSample[MAXPOPS];

        int i = 0;
        float sum = 0.0;
        int offset = ind*MAXPOPS;
        for(i = 0; i < MAXPOPS; i++){
            GammaSample[i] = RGammaDisc(Parameters[i],1,randState);
            sum += GammaSample[i];
        }
        for(i = 0; i < MAXPOPS; i++){
            TestQ[i+offset] = GammaSample[i]/sum;
        }
        saveRndDiscState(randState);
        ind += get_global_size(0);
    }
}


__kernel void FillArrayWRandom(
        __global float *randomArr,
        __global uint *randGens,
        const int length
        )
{
    int pos = get_global_id(0);
    mwc64x_state_t rng = getRandGen(randGens,pos);
    uint i;
    float val;
    ulong samplesPerstream = length/get_global_size(0);
    uint offset = pos*samplesPerstream;
    if (offset < length){
        for(i = 0; i < samplesPerstream && i+offset < length; i++){
            randomArr[offset+i] = uintToUnit(MWC64X_NextUint(&rng));
        }
        saveRandGen(randGens,pos,rng);
    }
}

__kernel void PopNormals(
        __global float *Prev,
        __global float *norms,
        __global uint *randGens,
        const float SD,
        const int length)
{
    RndDiscState randState[1];
    initRndDiscState(randState,randGens,0);
    int i;
    for( i=0; i < length; i+= 2){
        float2 rnorms = BoxMuller(randState);
        float oldf;
        oldf = Prev[i];
        norms[i] = rnorms.x*SD + oldf;

        if (i +1 < length){
            oldf = Prev[i+1];
            norms[i+1] = rnorms.y*SD + oldf;
        }
    }
    saveRndDiscState(randState);
}
