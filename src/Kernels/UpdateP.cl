__kernel void GetNumFromPops (
    __global int* Z, /* input */
    __global int* Geno,/* input */
    __global int* NumAlleles,/* input */
    __global int* popflags,/* input */
    __global int* NumAFromPops,/* output */
    __global int* error
)
{
    int ind = get_global_id(0);
    int loc = get_global_id(1);
    int numalleles = NumAlleles[loc];
    int offset = loc*MAXPOPS*MAXALLELES;
    int pos,line,popvalue,allelevalue;

    if(ind < NUMINDS && loc < NUMLOCI){
        if (!PFROMPOPFLAGONLY || popflags[ind] == 1) {    /*individual must have popflag turned on*/
            for (line = 0; line < LINES; line++) {
                popvalue = Z[ZPos (ind, line, loc)];
                allelevalue = Geno[GenPos (ind, line, loc)];

                if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
                    pos = NumAFromPopPos (popvalue, allelevalue) + offset;
                    atomic_add(&NumAFromPops[pos],1);
                }

            }
        }
    }
}

__kernel void UpdateP (
    __global double* P,
    __global double* LogP,
    __global double* lambda,
    __global int* NumAlleles,
    __global int* NumAFromPops,
    __global double* randomArr,
    __global int* error
)
{
    int loc = get_global_id(0);
    int numalleles = NumAlleles[loc];
    double Parameters[MAXALLELES];
    RndDiscState randState[1];
    int pop, allele;

    if(loc < NUMLOCI){
        int offset = loc*MAXPOPS*MAXALLELES;
        int pos,line,popvalue,allelevalue;
        initRndDiscState(randState,randomArr,MAXPOPS*MAXALLELES*MAXRANDOM);
        rndDiscStateReset(randState,offset);

        for (pop = 0; pop < MAXPOPS; pop++) {
            for (allele = 0; allele < numalleles; allele++) {
                    Parameters[allele] = lambda[pop] + NumAFromPops[NumAFromPopPos (pop, allele)+offset];
            }
            /*return a value of P simulated from the posterior Di(Parameters) */
            /*O(NumAlleles[loc]) */
            LogRDirichletDisc (Parameters, numalleles,
                               P + PPos (loc, pop, 0),
                               LogP +PPos(loc,pop,0),
                               randState);

        }
        if (randState->randomValsTaken > randState->maxrandom) {
            error[0] = KERNEL_OUT_OF_BOUNDS;
            error[1] = (randState->randomValsTaken - randState->maxrandom)/(MAXALLELES*MAXPOPS);
        }
    }
}
