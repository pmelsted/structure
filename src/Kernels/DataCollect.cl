void UpdateSumsPop (
        __global double *lambda,
        __global double *sumlambda,
        __global double *Fst,
        __global double *FstSum){
    int pop = get_global_id(0);
    sumlambda[pop] += lambda[pop];
    if (FREQSCORR) {
        FstSum[pop] += Fst[pop];
    }
}

/* Data collection paralell over pop */
__kernel void DataCollectPop(
       __global double *Alpha,
       __global double *AlphaSum,
       __global double *lambda,
       __global double *lambdaSum,
       __global double *Fst,
       __global double *FstSum
       )
{
    int pop = get_global_id(0);
    UpdateSumsPop(lambda,lambdaSum,Fst,FstSum);
    int pos,loc;
    if (LOCPRIOR && NOADMIX==0) {
            for (loc=0; loc<=NUMLOCATIONS; loc++) {
                pos = AlphaPos(loc, pop);
                AlphaSum[pos] += Alpha[pos];
            }
    } else if (!(NOADMIX) && (!(NOALPHA))) {
            AlphaSum[pop] += Alpha[pop];
    }

    /*if (pop == 0){
        if (LOCPRIOR)
            for (i=0; i<LocPriorLen; i++) {
                sumLocPrior[i] += LocPrior[i];
            }
    }*/
}

/* Data collection paralell over ind */
void UpdateSumsInd (
        __global double *Q,
        __global double *QSum,
        __global int *AncestDist
        )
{
    int pop = get_global_id(0);
    int ind = get_global_id(1);

    QSum[QPos (ind, pop)] += Q[QPos (ind, pop)];
    if (ANCESTDIST) {
        int box = ((int) (Q[QPos (ind, pop)] * ((double) NUMBOXES)));
        /*printf("%1.3f__%d  ",Q[QPos(ind,pop)],box); */
        if (box == NUMBOXES) {
            box = NUMBOXES - 1;    /*ie, Q = 1.000 */
        }
        AncestDist[AncestDistPos (ind, pop, box)]++;
    }
}

__kernel void DataCollectInd(
        __global double *Q,
        __global double *QSum,
        __global int *AncestDist
        )
{
    int pop = get_global_id(0);
    int ind = get_global_id(1);
    UpdateSumsInd(Q,QSum,AncestDist);
}

/* Data collection parallell over loc */
void UpdateSumsLoc (
        __global int *NumAlleles,
        __global double *P,
        __global double *PSum,
        __global double *Epsilon,
        __global double *SumEpsilon
        )
{
    int pop = get_global_id(0);
    int loc = get_global_id(1);
    int allele;

    for (allele = 0; allele < NumAlleles[loc]; allele++) {
        PSum[PPos (loc, pop, allele)] += P[PPos (loc, pop, allele)];
    }
    if (FREQSCORR){
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
            SumEpsilon[EpsPos (loc, allele)] += Epsilon[EpsPos (loc, allele)];
        }
    }
}

__kernel void DataCollectLoc(
        __global int *NumAlleles,
        __global double *P,
        __global double *PSum,
        __global double *Epsilon,
        __global double *SumEpsilon
        ) 
{
    int pop = get_global_id(0);
    int loc = get_global_id(1);
    UpdateSumsLoc(NumAlleles,P,PSum,Epsilon,SumEpsilon);
}

__kernel void ComputeProbFinish(
        __global double *loglike,
        __global double *sumlikes,
        __global double *sumsqlikes
        )
{
    double like = loglike[0];
    sumlikes[0] += like;
    sumsqlikes[0] += like*like;
}
