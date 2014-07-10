#include "Kernels/randGen.cl"

__kernel void Dirichlet(
        __global double *Parameters,
        __global uint *randGens,
        __global double *TestQ)
{
    int ind = get_global_id(0);
    RndDiscState randState[1];

    while (ind < NUMINDS){
        initRndDiscState(randState,randGens,ind);
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
        saveRndDiscState(randState);
        ind += get_global_size(0);
    }
}


__kernel void FillArrayWRandom(
        __global double *randomArr,
        __global uint *randGens,
        const int length
        )
{
    int pos = get_global_id(0);
    mwc64x_state_t rng = getRandGen(randGens,pos);
    uint i;
    double val;
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
        __global double *Prev,
        __global double *norms,
        __global uint *randGens,
        const double SD)
{
    int pop = get_global_id(0);
    while(pop < MAXPOPS){
        RndDiscState randState[1];
        initRndDiscState(randState,randGens,pop);
        double2 rnorms = BoxMuller(randState);
        double oldf = Prev[pop];
        norms[pop] = rnorms.x*SD + oldf;
        saveRndDiscState(randState);
        pop += get_global_size(0);
    }
}
