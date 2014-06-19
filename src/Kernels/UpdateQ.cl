__kernel void UpdateQMetro (
    __global int *Geno,
    __global int *PreGeno,
    __global double *Q,
    __global double *P,
    __global double *Alpha,
    __global int *individualLocs,
    __global int *individualPopFlags,
    __global int *Recessive,
    __global double * randomArr,
    const int rep,
    __global int* error
)
{
    double CurrentQ[MAXPOPS];             /*[MAXPOPS]; */
    double TestQ[MAXPOPS];                /*[MAXPOPS]; */
    double PriorQ1[MAXPOPS];
    int pop;
    double logdiff;
    double randomnum;
    int apos;
    int numhits = 0;


    RndDiscState randState[1];
    int ind = get_global_id(0);

    if (ind < NUMINDS) {
        if (!((USEPOPINFO) && (individualPopFlags[ind]))) {
            rndDiscStateReset(randState, ind*MAXRANDOM);
        }
        #if LOCPRIOR
        apos = AlphaPos(individualLocs[ind],0);
        for (pop = 0; pop < MAXPOPS; pop++) {
            PriorQ1[pop] = Alpha[apos+pop];
        }
        #else
        for (pop = 0; pop < MAXPOPS; pop++) {
            PriorQ1[pop] = Alpha[pop];
        }
        #endif

        RDirichletDisc (PriorQ1, MAXPOPS, TestQ,
                        randState);     /*return TestQ, sampled from the prior */



        if (rep == 0) {
            for (pop = 0; pop < MAXPOPS; pop++) {
                Q[QPos (ind, pop)] = (double) 1 / MAXPOPS;
            }
        }

        for (pop = 0; pop < MAXPOPS; pop++) {
            CurrentQ[pop] = Q[QPos (ind, pop)];
        }


        logdiff = 0.0;
        /* logdiff += log(TestQ[pop]) - log(CurrentQ[pop]);
        logdiff = logdiff*(alpha-1.0); removed prior prob bit */

        logdiff += CalcLikeInd (Geno, PreGeno, TestQ, P, ind,
                                Recessive);  /*likelihood bit */
        logdiff -= CalcLikeInd (Geno, PreGeno, CurrentQ, P, ind,
                                Recessive);

        randomnum = rndDisc(randState);
        if (randomnum < exp (logdiff)) {    /*accept */
            for (pop = 0; pop < MAXPOPS; pop++) {
                Q[QPos (ind, pop)] = TestQ[pop];
            }
            numhits++;
        }
    }
}
