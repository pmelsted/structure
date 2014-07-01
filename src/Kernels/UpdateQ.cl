__kernel void MetroAcceptTest(
        __global double *TestQ,
        __global double *Q,
        __global double *randomArr,
        __global double *logdiffs,
        __global int *popflags)
{
    int ind = get_global_id(0);
    int pop;
    RndDiscState randState[1];

    initRndDiscState(randState,randomArr,1);
    rndDiscStateReset(randState,NUMINDS*MAXRANDOM+ind);
    if (!((USEPOPINFO) && (popflags[ind]))) {
        if(rndDisc(randState) < exp(logdiffs[ind])){
            for (pop = 0; pop < MAXPOPS; pop++) {
                Q[QPos (ind, pop)] = TestQ[QPos(ind,pop)];
            }
        }
    }
}

__kernel void GetNumLociPops(
        __global int *Z,
        __global int *popflags,
        __global int *NumLociPops)
{
    int ind = get_global_id(0);
    int loc = get_global_id(1);
    int offset = ind*MAXPOPS;
    int line, from,pop;
    /* initialize the NumLociPops array */
    if(loc == 0){
        for(pop = 0; pop < MAXPOPS; pop++){
            NumLociPops[pop+offset] = 0;
        }
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
    if(ind < NUMINDS && loc < NUMLOCI) {
        if (!((USEPOPINFO) && (popflags[ind]))) {
            for (line = 0; line < LINES; line++) {
                from = Z[ZPos (ind, line, loc)];
                if (from != UNASSIGNED) {
                    atomic_add(&NumLociPops[from+offset],1);
                }
            }
        }
    }
}

__kernel void UpdQDirichlet(
        __global double *Alpha,
        __global int *NumLociPops,
        __global double *randomArr,
        __global double *Q)
{
    int ind = get_global_id(0);
    RndDiscState randState[1];

    initRndDiscState(randState,randomArr,MAXRANDOM);
    rndDiscStateReset(randState,ind*MAXRANDOM);
    double GammaSample[MAXPOPS];

    int i = 0;
    double sum = 0.0;
    int offset = ind*MAXPOPS;
    double param;
    for(i = 0; i < MAXPOPS; i++){
        param = Alpha[i]+NumLociPops[i+offset];
        GammaSample[i] = RGammaDisc(param,1,randState);
        sum += GammaSample[i];
    }
    for(i = 0; i < MAXPOPS; i++){
        Q[i+offset] = GammaSample[i]/sum;
    }
}


