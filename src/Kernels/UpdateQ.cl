/*__kernel void UpdateQMetro (
    __global int *Geno,
    __global int *PreGeno,
    __global double *Q,
    __global double *P,
    __global double *Alpha,
    [>__global struct IND *Individual,<]
    __global int *individualLocs,
    __global int *individualPopFlags,
    __global int *Recessive,
    __global double * randomArr,
    __global int* error
    )
{
    double CurrentQ[MAXPOPS];             [>[MAXPOPS]; <]
    double TestQ[MAXPOPS];                [>[MAXPOPS]; <]
    double PriorQ1[MAXPOPS];
    int pop;
    double logdiff;
    double randomnum;
    int apos;
    int numhits = 0;


    RndDiscState randState[1];
    int ind = get_global_id(0);

    if (ind < NUMINDS){
        if (!((USEPOPINFO) && (individualPopFlags[ind]))) {
            rndDiscStateReset(randState, ind*MAXRANDOM);
        }
        #if LOCPRIOR
            apos = AlphaPos(individualLocs[ind],0);
            for (pop = 0; pop < MAXPOPS; pop++){
                PriorQ1[pop] = Alpha[apos+pop];
            }
        #else
            for (pop = 0; pop < MAXPOPS; pop++){
                PriorQ1[pop] = Alpha[pop];
            }
        #endif

      RDirichletDisc (PriorQ1, MAXPOPS, TestQ,randState);     [>return TestQ, sampled from the prior <]
    }
}*/
