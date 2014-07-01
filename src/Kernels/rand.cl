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

