#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"
#include "../output.h"
#include "../Kernels.h"
#include "../CalcLike.h"
#include "ForwardAndBackward.h"

void UpdateQAdmixture (float *Q, int *Z, float *Alpha,
                       struct IND *Individual, float * randomArr);


/*------------------------------------------*/
void GetNumLociPop (int *NumLociPop, int *Z, int ind)
{
    /*Fill in the number of alleles that each individual has from each pop */
    int loc, line, pop;
    int from;

    for (pop = 0; pop < MAXPOPS; pop++) {
        NumLociPop[pop] = 0;
    }

    for (loc = 0; loc < NUMLOCI; loc++) {
        for (line = 0; line < LINES; line++) {
            from = Z[ZPos (ind, line, loc)];
            if (from != UNASSIGNED) {
                NumLociPop[from]++;
            }
        }
    }
}

float mapLogDiffsFunc(float *Q,float *TestQ,
                       float *P,int *Geno, int ind, int loc){
    int allele, line, pop;
    float termP;
    float termM;
    float logterm = 0.0;

    for (line = 0; line < LINES; line++) {
        allele = Geno[GenPos (ind, line, loc)];
        if (allele != MISSING) {
            termP = 0.0;
            termM = 0.0;
            for (pop = 0; pop < MAXPOPS; pop++) {
                termP += TestQ[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                termM += Q[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
            }
            logterm += log(termP) - log(termM);
        }
    }
    return logterm;
}

float verifyReduction(float *Q,float *TestQ,
                       float *P,int *Geno, float * logdiffs){

    float * cpulogdiffs;
    int ind, loc;
    float logdiff;
    float diff;
    float ret = 0.0;
    cpulogdiffs = calloc(NUMINDS,sizeof(float));

    for(ind = 0; ind < NUMINDS; ind++){
        logdiff = 0.0;
        for(loc = 0; loc < NUMLOCI; loc++){
            logdiff += mapLogDiffsFunc(Q,TestQ,P,Geno,ind,loc);
        }
        cpulogdiffs[ind] = logdiff;
    }

    for(ind = 0; ind < NUMINDS; ind++){
        diff = fabs(logdiffs[ind] -cpulogdiffs[ind]);
        if ( diff > 10e-6){
            printf("Reduction error: ind: %d diff: is %f!\n",ind,diff);
            printf("cl %.10f, cpu %.10f\n",logdiffs[ind],cpulogdiffs[ind]);
            ret = diff;
        }
    }

    return ret;

}

/*----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate locprior*/
void UpdateQMetroCL (CLDict *clDict,int *Geno, int *PreGeno, float *Q, float *P,
                   float *Alpha, int rep, struct IND *Individual, int *Recessive,
                   float * randomArr)

{
    /*int pop;*/
    /*float randomnum;*/
    /*int ind;*/
    /*int numhits = 0;*/

    /*float *TestQ;*/
    /*float *logdiffs;*/
    /*float *logterms;*/

    /*RndDiscState randState[1];*/
    size_t global[2];

    /*
     * Removed:
     *   RECESSIVEALLELES
     *   LOCPRIOR
     *
     */
    /*float verif;*/
    /*logterms = calloc(NUMINDS*NUMLOCI,sizeof(float));*/
    /*TestQ = calloc (NUMINDS*MAXPOPS, sizeof (float));*/
    /*logdiffs = calloc(NUMINDS,sizeof(float));*/


    /*initRndDiscState(randState,randomArr,NUMINDS + NUMINDS*MAXRANDOM);*/

    /*writeBuffer(clDict,randomArr,sizeof(float) * (NUMINDS*MAXRANDOM+NUMINDS),RANDCL,"Random");*/

    /* ======== Sample ====== */
    global[0] = NUMINDS;
    runKernel(clDict,RDirichletSampleKernel,1,global,"Dirichlet");

    /*readBuffer(clDict,TestQ,sizeof(float) *QSIZE,TESTQCL,"TestQ");*/

    /*for (ind = 0; ind < NUMINDS; ind++) {*/
        /*RDirichletDisc(Alpha, MAXPOPS, TestQ,ind*MAXPOPS,randState);*/
    /*}*/



    /* ======== Calculate likelihood ====== */
    global[0] = NUMLOCI;
    global[1] = NUMINDS;

    runKernel(clDict,mapReduceLogDiffsKernel,2,global,"reduceLogDiffs");

    /*
    logterms = calloc(NUMINDS*NUMLOCI,sizeof(float));
    TestQ = calloc (NUMINDS*MAXPOPS, sizeof (float));
    logdiffs = calloc(NUMINDS,sizeof(float));
    readBuffer(clDict,TestQ,sizeof(float) *QSIZE,TESTQCL,"TestQ");
    readBuffer(clDict,logdiffs,sizeof(float) * NUMINDS,LOGDIFFSCL,"Logdiffs");
    readBuffer(clDict,Q,sizeof(float) *QSIZE,QCL,"Q");
    readBuffer(clDict,P,sizeof(float) *PSIZE,PCL,"P");
    readBuffer(clDict,Geno,sizeof(int) *GENOSIZE,GENOCL,"Geno");

    verif  = verifyReduction(Q,TestQ,P,Geno,logdiffs);

    if (verif > 10e-6){
        handleCLErr(1,clDict,"red err");
    }
    free (TestQ);
    free (logdiffs);
    free(logterms);
    */

    /*CalcLogdiffsCL(clDict,Geno,TestQ,Q,P,logdiffs);*/

    /* ========= Acceptance test ========= */
    global[0] = NUMINDS;
    runKernel(clDict,MetroAcceptTestKernel,1,global,"MetroAcceptTest");


    /*for (ind = 0; ind < NUMINDS; ind++){
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            randomnum = rndDisc(randState);
            if (randomnum < exp (logdiffs[ind])) {    [>accept <]
                for (pop = 0; pop < MAXPOPS; pop++) {
                    Q[QPos (ind, pop)] = TestQ[QPos(ind,pop)];
                }
                numhits++;
            }
        }
    }*/

    /* TODO: count hits in Acceptance test (put 0,1 in array and reduce?)*/
    /*if (REPORTHITRATE) {           [>does this once every UPDATEFREQ reps <]*/
        /*if ((int) (rep - METROFREQ) / UPDATEFREQ < (int) rep / UPDATEFREQ) {*/
            /*printf ("Acceptance rate in UpdateQMetro %1.3f\n",*/
                    /*(float) numhits / NUMINDS);*/
        /*}*/
    /*}*/

    /*free (TestQ);*/
    /*free (logdiffs);*/
    /*free(logterms);*/
}



/*----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate locprior*/
void UpdateQMetro (int *Geno, int *PreGeno, float *Q, float *P,
                   float *Alpha, int rep, struct IND *Individual, int *Recessive,
                   float * randomArr)
/*
 * The goal with this function is to improve mixing when alpha is
 * small.  The problem is that in that situation the Q's can't move
 * efficiently.  What I do here is simply pick a proposal Q from the
 * prior, and try to move there.  I move from q0->q1 with prob: Min(1,
 * [P(X|q1,P)/P(X|q0,P)]
 * --- Notes 28 July 99
 */

/*
 * note that this code overlaps substantially with updateqmetrorecombine (the
 * version used for the linkage model).  Both versions should be updated in tandem
 */

{
    float *PriorQ1;    /*[MAXPOPS]; */
    float *CurrentQ;             /*[MAXPOPS]; */
    float *TestQ;                /*[MAXPOPS]; */
    int pop;
    float logdiff;
    float randomnum;
    int ind;
    int numhits = 0;

    RndDiscState randState[1];
    /*  int i, ok; */

    /*  PriorQ1 = calloc (MAXPOPS, sizeof (float)); */
    CurrentQ = calloc (MAXPOPS, sizeof (float));
    TestQ = calloc (MAXPOPS, sizeof (float));

    if ((CurrentQ == NULL) || (TestQ == NULL)) {
        printf ("WARNING: error in assigning memory in function UpdateQMetro\n");
        Kill ();
    }

    /*  for (pop = 0; pop < MAXPOPS; pop++)
        PriorQ1[pop] = Alpha[pop];
    */

    PriorQ1 = Alpha;

    initRndDiscState(randState,randomArr,MAXRANDOM);
    for (ind = 0; ind < NUMINDS; ind++) {
        rndDiscStateReset(randState, ind*MAXRANDOM);
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            /* ie don't use individuals for whom prior pop info is used */

            /*-------compute/record newq and oldq----------*/
            /*ok = 1;
            do
            { */

            if (LOCPRIOR) {
                PriorQ1 = &Alpha[AlphaPos(Individual[ind].myloc, 0)];
            }
            /*return TestQ, sampled from the prior */
            RDirichletDisc (PriorQ1, MAXPOPS, TestQ,0, randState);

            /*  for (i=0;i<MAXPOPS;i++)
            if (TestQ[i]==0) { ok=0; break;}
            }
            while (ok==0); */

            for (pop = 0; pop < MAXPOPS; pop++) {
                CurrentQ[pop] = Q[QPos (ind, pop)];
            }

            /*-------Do metropolis test of newq-------*/

            logdiff = 0.0;
            /* logdiff += log(TestQ[pop]) - log(CurrentQ[pop]);
            logdiff = logdiff*(alpha-1.0); removed prior prob bit */

            if (LINES==2 && RECESSIVEALLELES) {
                logdiff += CalcLikeInd (Geno, PreGeno, TestQ, P, ind,
                                        Recessive);  /*likelihood bit */
                logdiff -= CalcLikeInd (Geno, PreGeno, CurrentQ, P, ind, Recessive);
            } else {
                logdiff = CalcLikeIndDiff (Geno, PreGeno, TestQ, CurrentQ, P, ind, Recessive);
            }

            randomnum = rndDisc(randState);
            if (randomnum < exp (logdiff)) {    /*accept */
                for (pop = 0; pop < MAXPOPS; pop++) {
                    Q[QPos (ind, pop)] = TestQ[pop];
                }
                numhits++;
            }
            /*for (pop=0;pop<MAXPOPS; pop++)
            printf("%1.3f %1.3f    ",CurrentQ[pop],TestQ[pop]);
            if (randomnum < exp(logdiff) ) printf("  Accepted ");
            printf("\n"); */
        }
    }
    if (REPORTHITRATE) {           /*does this once every UPDATEFREQ reps */
        if ((int) (rep - METROFREQ) / UPDATEFREQ < (int) rep / UPDATEFREQ) {
            printf ("Acceptance rate in UpdateQMetro %1.3f\n",
                    (float) numhits / NUMINDS);
        }
    }

    /*  free (PriorQ1); */
    free (CurrentQ);
    free (TestQ);
}


/*----------------------------------------*/
/*O(NUMINDS*MAXPOPS*LINES*NUMLOCI)*/
void
UpdateQNoAdmix (int *Geno, float *Q, float *P, struct IND *Individual,
                float *LocPrior,float * randomArr)
/*
 * Assign each individual to exactly one population according to the
 * conditional probabilities.
 */
{

    int ind, line, loc, pop;
    int allele;
    float *ProbsVector;          /*[MAXPOPS] */
    float sumlogs, sum;
    float runningtotal;
    float max=0.0, prob;
    int pickedpop;

    ProbsVector = calloc (MAXPOPS, sizeof (float));

    for (ind = 0; ind < NUMINDS; ind++) {
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            /* ie don't use individuals for whom prior pop info is used */

            /*make a vector of log probs for each pop */
            for (pop = 0; pop < MAXPOPS; pop++) {
                sumlogs = 0;

                /*Melissa added 7/12/07*/
                if (LOCPRIOR) {
                    prob = LocPrior[LocPriorPos(Individual[ind].myloc, pop)];
                } else {
                    prob = 1.0;
                }

                runningtotal = prob;
                for (line = 0; line < LINES; line++) {
                    for (loc = 0; loc < NUMLOCI; loc++) {
                        allele = Geno[GenPos (ind, line, loc)];
                        if (allele != MISSING) {
                            runningtotal *= P[PPos (loc, pop, allele)];
                            if (runningtotal < UNDERFLO) {
                                /*this is to avoid having to
                                                take logs all the time */
                                sumlogs += log (runningtotal);
                                runningtotal = 1;
                            }
                        }
                    }
                }
                ProbsVector[pop] = sumlogs + log (runningtotal);
                if (pop==0 || ProbsVector[pop] > max) {
                    max = ProbsVector[pop];
                }
            }

            sum = 0.0;
            for (pop=0; pop < MAXPOPS; pop++) {
                sum += (ProbsVector[pop] = exp(ProbsVector[pop]-max));
            }

            for (pop = 0; pop < MAXPOPS; pop++) {
                Q[QPos (ind, pop)] = 0.0;
            }

            pickedpop = PickAnOption (MAXPOPS, sum, ProbsVector);
            Q[QPos (ind, pickedpop)] = 1.0;

        }
    }

    free (ProbsVector);
}


/*-----------------------------------------*/
void UpdateQAdmixture (float *Q, int *Z, float *Alpha,
                       struct IND *Individual,float * randomArr)
/*update Q: proportion of ancest of each ind in each pop. */
{
    int *NumLociPop;              /*[MAXPOPS] -- number of alleles from each pop (by ind) */
    float *Parameters;           /*[MAXPOPS] -- dirichlet parameters of posterior on Q */
    int ind, pop, loc;
    float *usealpha=NULL;

    NumLociPop = calloc (MAXPOPS, sizeof (int));
    Parameters = calloc (MAXPOPS, sizeof (float));
    if ((NumLociPop == NULL) || (Parameters == NULL)) {
        printf ("WARNING: unable to allocate array space in UpdateQAdmixture\n");
        Kill ();
    }

    if (LOCPRIOR==0) {
        usealpha=Alpha;
    }
    for (ind = 0; ind < NUMINDS; ind++) {
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            /* ie don't use individuals for whom prior pop info is used */
            GetNumLociPop (NumLociPop, Z, ind);
            if (LOCPRIOR) {
                loc = Individual[ind].myloc;
                usealpha = &Alpha[AlphaPos(loc, 0)];
            }

            for (pop = 0; pop < MAXPOPS; pop++) {       /*compute dirichlet parameters */
                Parameters[pop] = usealpha[pop] + NumLociPop[pop];
            }
            /*generate dirichlet RVs */
            RDirichlet (Parameters, MAXPOPS, Q + QPos (ind, 0));
        }
    }
    free (NumLociPop);
    free (Parameters);
}



/*-----------------------------------------*/
void UpdateQAdmixtureCL (CLDict *clDict,float *Q, int *Z, float *Alpha,
                       struct IND *Individual,float * randomArr)
/*update Q: proportion of ancest of each ind in each pop. */
{
    /*int *NumLociPop;              [>[MAXPOPS] -- number of alleles from each pop (by ind) <]
    float *Parameters;           [>[MAXPOPS] -- dirichlet parameters of posterior on Q <]
    int ind, pop, loc;
    float *usealpha=NULL;
    int offset;*/

    /*int *NumLociPops;*/
    size_t global[2];
    /*NumLociPop = calloc (MAXPOPS, sizeof (int));*/
    /*Parameters = calloc (MAXPOPS, sizeof (float));*/

    /*NumLociPops = calloc (NUMINDS*MAXPOPS, sizeof (int));*/
    /*if ((NumLociPop == NULL) || (Parameters == NULL)) {
        printf ("WARNING: unable to allocate array space in UpdateQAdmixture\n");
        Kill ();
    }*/

    /*if (LOCPRIOR==0) {
        usealpha=Alpha;
    }*/
    /* Clear the buffer */
    /*writeBuffer(clDict,NumLociPops,sizeof(int)* NUMINDS*MAXPOPS, NUMLOCIPOPSCL,"NumLociPops");*/
    global[0] = NUMINDS;
    global[1] = NUMLOCI;

    runKernel(clDict,GetNumLociPopsKernel,2,global,"getNumLociPops");
    /*readBuffer(clDict,NumLociPops,sizeof(int)* NUMINDS*MAXPOPS, NUMLOCIPOPSCL,"NumLociPops");*/

    /*writeBuffer(clDict,randomArr,sizeof(float) *NUMINDS*MAXRANDOM,RANDCL,"Random");*/
    runKernel(clDict,UpdQDirichletKernel,1,global,"UpdateQDirichlet");
    /*GetNumLociPopsCL (clDict,NumLociPops, Z, ind);*/
    /*for (ind = 0; ind < NUMINDS; ind++) {
        offset = ind*MAXPOPS;
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            [> ie don't use individuals for whom prior pop info is used <]
            [>GetNumLociPop (NumLociPop, Z, ind);<]
            [>for (pop = 0; pop < MAXPOPS; pop++) {<]
                [>if(NumLociPop[pop] != NumLociPops[pop+offset]){<]
                    [>ReleaseCLDict(clDict);<]
                    [>printf("%d, %d\n",NumLociPop[pop],NumLociPops[pop+offset]);<]
                    [>exit(EXIT_FAILURE);<]
                [>}<]
            [>}<]
            if (LOCPRIOR) {
                loc = Individual[ind].myloc;
                usealpha = &Alpha[AlphaPos(loc, 0)];
            }

            for (pop = 0; pop < MAXPOPS; pop++) {       [>compute dirichlet parameters <]
                Parameters[pop] = usealpha[pop] + NumLociPops[pop+offset];
            }
            [>generate dirichlet RVs <]
            RDirichlet (Parameters, MAXPOPS, Q + QPos (ind, 0));
        }
    }
    free (NumLociPop);
    free (Parameters);*/


    /*free(NumLociPops);*/
}



/*-----------------------------------------*/
void
UpdateQWithPops (int *Geno, float *Q, float *P, int *Z, float *Alpha,
                 int rep, struct IND *Individual, float *UsePopProbs,
                 float * randomArr)
/*
 * This version of UpdateQ is used when the USEPOPINFO flag is
 * turned on, indicating that the prior information about population
 * information should be taken into account.  It assumes that individuals
 * are normally from the given populations, but are occasionally
 * migrants, or children, grandchildren, etc of migrants.  Individuals
 * can also be marked as unknown, in which case the usual method for
 * estimating Q is used. Individuals who lack prior pop info are treated
 * by the regular Q updates
 */
{
    float *Input, *Output;       /*both [MAXPOPS] */
    int ind, pop, loc, allele1, allele2, homepop;
    int pop1;
    int gen;
    float sumlogs;
    float sumsofar;
    float runningtotal;
    float *PostProbs;            /*[MAXPOPS][GENSBACK+1] */
    float *Priors;               /*[GENSBACK+1] */
    float p;
    float weights;
    float sumweights;
    float sumprobs;
    float rannum;
    float sqrtunder = sqrt (UNDERFLO);
    float like=0.0;              /*likelihood of current genotype */
    float maxprob;


    /*compute prior probs on each option */
    /*the correct-as-assigned option has weight MIGRPRIOR; the remaining
      probability is split among the other options */

    Input = calloc (MAXPOPS, sizeof (float));
    Output = calloc (MAXPOPS, sizeof (float));
    PostProbs = calloc (MAXPOPS * (GENSBACK + 1), sizeof (float));
    Priors = calloc (GENSBACK + 1, sizeof (float));
    if ((Input == NULL) || (Output == NULL) || (PostProbs == NULL)
            || (Priors == NULL)) {
        printf ("Warning!! Error in assigning memory in UpdateQWithPops\n");
        Kill ();
    }

    weights = 0.5;
    sumweights = 0.0;
    for (gen = 0; gen < GENSBACK + 1; gen++) {
        weights *= 2;               /*sequence goes 1, 2, 4, 8,.... */
        Priors[gen] = weights;
        sumweights += weights * (MAXPOPS - 1);
    }

    for (gen = 0; gen < GENSBACK + 1; gen++) {
        Priors[gen] = MIGRPRIOR * Priors[gen] / sumweights;
    }


    /*now compute ancestry */
    for (ind = 0; ind < NUMINDS; ind++) {
        if (Individual[ind].PopFlag) {        /*use prior population information */
            homepop = Individual[ind].Population - 1;
            for (gen = 0; gen < GENSBACK + 1; gen++) {
                if (gen>0) {
                    p = exp ((gen - 1) * log (0.5));
                } else {
                    p=1.0;
                }

                for (pop = 0; pop < MAXPOPS; pop++) {
                    sumlogs = 0;
                    runningtotal = 1;

                    for (loc = 0; loc < NUMLOCI; loc++) {
                        allele1 = Geno[GenPos (ind, 0, loc)];
                        if (LINES>1) {
                            allele2 = Geno[GenPos (ind, 1, loc)];
                        } else {
                            allele2=MISSING;
                        }

                        if ((allele1 != MISSING) && (allele2 != MISSING)) {
                            /* no missing alleles */
                            if (gen == 0) {             /*both alleles from same population */
                                like = P[PPos (loc, pop, allele1)] * P[PPos (loc, pop, allele2)];
                            } else if (gen > 0) {
                                /*at most one allele from the outside population */
                                like =
                                    0.5 * p * P[PPos (loc, pop, allele1)] * P[PPos (loc, homepop, allele2)] +
                                    0.5 * p * P[PPos (loc, homepop, allele1)] * P[PPos (loc, pop, allele2)] +
                                    (1.0 - p) * P[PPos (loc, homepop, allele1)] * P[PPos (loc, homepop, allele2)];
                            }
                        } else if (allele1 != MISSING) {       /*1 allele missing */
                            like = p * P[PPos (loc, pop, allele1)]
                                   + ((float) 1.0 - p) * P[PPos (loc, homepop, allele1)];
                        } else if (allele2 != MISSING) {
                            like = p * P[PPos (loc, homepop, allele2)]
                                   + ((float) 1.0 - p) * P[PPos (loc, homepop, allele2)];
                        } else {
                            like = 1.0;       /*both alleles missing */
                        }

                        if (like > sqrtunder) {      /*0-values lead to underflow */
                            runningtotal *= like;
                        } else {
                            runningtotal *= sqrtunder;
                        }

                        if (runningtotal < sqrtunder) {
                            /*this is to avoid having to
                                                                  take logs all the time */
                            sumlogs += log (runningtotal);
                            runningtotal = 1;
                        }
                    }

                    /*note: this next bit of code is changed from
                      previously, where I took off some standard
                      normalizing constant, and put a real probability
                      (not log prob) into PostProbs */
                    sumlogs = sumlogs + log (runningtotal);
                    PostProbs[PostProbsPos (pop, gen)] = sumlogs;

                }
            }

            /*find maximum log prob */
            maxprob = PostProbs[PostProbsPos (0, 0)];
            for (gen = 0; gen < GENSBACK + 1; gen++) {
                for (pop = 0; pop < MAXPOPS; pop++) {
                    if (PostProbs[PostProbsPos (pop, gen)] > maxprob) {
                        maxprob = PostProbs[PostProbsPos (pop, gen)];
                    }
                }
            }

            /*subtract off minimum to avoid underflow-type problems,and
              exponentiate to turn into regular probabilities */
            for (gen = 0; gen < GENSBACK + 1; gen++) {
                for (pop = 0; pop < MAXPOPS; pop++) {
                    PostProbs[PostProbsPos (pop, gen)] =
                        exp (PostProbs[PostProbsPos (pop, gen)] - maxprob);
                }
            }
            /*compute posterior probabilities (unnormalized) */
            sumprobs = 0.0;
            for (gen = 0; gen < GENSBACK + 1; gen++) {
                for (pop = 0; pop < MAXPOPS; pop++) {
                    if ((pop == homepop) && (gen == 0)) {   /* non-migrant case */
                        PostProbs[PostProbsPos (pop, gen)] *= (1 - MIGRPRIOR);
                        sumprobs += PostProbs[PostProbsPos (pop, gen)];
                    } else if ((pop == homepop) && (gen != 0)) {     /* dealt with elsewhere */
                        PostProbs[PostProbsPos (pop, gen)] = 0.0;
                    } else {
                        /*migrant case */
                        PostProbs[PostProbsPos (pop, gen)] =
                            PostProbs[PostProbsPos (pop, gen)] * Priors[gen];
                        sumprobs += PostProbs[PostProbsPos (pop, gen)];
                    }
                    /*printf("%2.8f  ",Infer1Probs[pop][gen]); */
                }
                /*printf("\n"); */
            }
            /*printf("%d %f\n\n",homepop+1,sumprobs); */

            if (rep >= BURNIN) {
                for (gen = 0; gen < GENSBACK + 1; gen++) {
                    for (pop = 0; pop < MAXPOPS; pop++) {
                        UsePopProbs[UsePPrPos (ind, pop, gen)]
                        += PostProbs[PostProbsPos (pop, gen)] / sumprobs;
                    }
                }
            }

            /*now pick an option from the posterior... */
            rannum = RandomReal (0.0, sumprobs);
            sumsofar = 0.0;
            for (gen = 0; gen < GENSBACK + 1; gen++) {
                for (pop = 0; pop < MAXPOPS; pop++) {
                    sumsofar += PostProbs[PostProbsPos (pop, gen)];
                    if (rannum < sumsofar) {        /*set ancest to this guy */
                        for (pop1 = 0; pop1 < MAXPOPS; pop1++) {
                            Q[QPos (ind, pop1)] = 0;
                        }
                        p = exp (gen * log (0.5));  /*amount of migrant ancestry */
                        Q[QPos (ind, pop)] = p;
                        Q[QPos (ind, homepop)] += 1 - p;
                        /*printf("P: %d, G: %d\n",pop+1,gen);
                          printf("%1.3f %1.3f %1.3f\n",
                          Ancest[ind][0],Ancest[ind][1],Ancest[ind][2]); */

                        /*set remaining ancestry to home pop */
                        gen = GENSBACK + 1;
                        pop = MAXPOPS;      /*breaking out of the loop */
                    }
                }
                /*printf("\n"); */
            }
        }
    }

    free (Input);
    free (Output);
    free (PostProbs);
    free (Priors);
}



/*-----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate LocPriors*/
void UpdateQ (int *Geno, int *PreGeno, float *Q, float *P, int *Z,
              float *Alpha, int rep,
              struct IND *Individual, float *UsePopProbs, int *Recessive, float *LocPrior,
              float * randomArr)
/*
 * update Q: proportion of ancest of each ind in each pop. Three
 * main options here.  One is update Q with admixture, one is update
 * Q with no admixture, and one is to use prior population info.
 * Even when USEPOPINFO is turned on, the other Q updates are still
 * called, in order to deal with any individuals who lack population
 * info.  The other functions are written to ignore individuals with
 * pop data when USEPOPINFO is on.
 */
{

    if (USEPOPINFO) {               /*update with prior population information */
        UpdateQWithPops (Geno, Q, P, Z, Alpha, rep, Individual, UsePopProbs,randomArr);
    }

    if (NOADMIX) {                  /*no admixture model */
        /* don't use ADMIBURNIN with LOCPRIOR models */
        if ((rep > ADMBURNIN) || (rep > BURNIN) || LOCPRIOR) {
            UpdateQNoAdmix (Geno, Q, P, Individual, LocPrior,randomArr);
        } else {
            UpdateQAdmixture (Q, Z, Alpha, Individual,randomArr);       /*initial burnin */
        }
    } else {
        /*admixture model */
        if (METROFREQ > 0 && rep%METROFREQ==0) {
            UpdateQMetro (Geno, PreGeno, Q, P, Alpha, rep, Individual, Recessive,
                          randomArr);
        } else {
            UpdateQAdmixture (Q, Z, Alpha, Individual,randomArr);
        }
    }
}

void
UpdateQMetroRecombine (int *Geno, float *Q, int *Z, float *P,
                       float *Alpha, int rep, struct IND *Individual,
                       float *Mapdistance, float *R,
                       float *Phase,int *Phasemodel,
                       float * randomArr)
/*
 * This function does the same job as UpdateQMetro in the case
 * when there is recombination.
 * the code is very similar, but a new routine is nevertheless a good idea
 */
{
    float *PriorQ;               /*[MAXPOPS]; */
    float *CurrentQ;             /*[MAXPOPS]; */
    float *TestQ;                /*[MAXPOPS]; */
    int pop;
    float logdiff;
    float randomnum;
    int ind;
    int numhits = 0;
    /*  int i, ok; */
    float *RTransitProb;
    float *IndividualQ;

    if (PHASED) {
        RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (float));
    } else {
        RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (float));
    }

    PriorQ = calloc (MAXPOPS, sizeof (float));
    CurrentQ = calloc (MAXPOPS, sizeof (float));
    TestQ = calloc (MAXPOPS, sizeof (float));
    IndividualQ = calloc (MAXPOPS, sizeof (float));

    if ((PriorQ == NULL) || (CurrentQ == NULL) || (TestQ == NULL)) {
        printf ("WARNING: error in assigning memory in function UpdateQMetro\n");
        Kill ();
    }

    for (pop = 0; pop < MAXPOPS; pop++) {
        PriorQ[pop] = Alpha[pop];
    }

    for (ind = 0; ind < NUMINDS; ind++) {
        if ((USEPOPINFO)
                && (Individual[ind].PopFlag)) {/*set Q for inds with prior info*/
            for (pop=0; pop<MAXPOPS; pop++) {
                if (pop==Individual[ind].Population-1) {
                    Q[QPos(ind,pop)]=1.0;
                } else {
                    Q[QPos(ind,pop)]=0.0;
                }
            }
        } else {
            /* ie don't use individuals for whom prior pop info is used */

            /*-------compute/record newq and oldq----------*/
            /*ok = 1;
              do
              { */
            RDirichlet (PriorQ, MAXPOPS,
                        TestQ);      /*return TestQ, sampled from the prior */
            /*  for (i=0;i<MAXPOPS;i++)
                if (TestQ[i]==0) { ok=0; break;}
                }
                while (ok==0); */


            if (rep ==
                    0) {             /*If this is the first rep, Q will not be initialized */
                for (pop = 0; pop < MAXPOPS; pop++) {
                    Q[QPos (ind, pop)] = (float) 1 / MAXPOPS;
                }
            }

            for (pop = 0; pop < MAXPOPS; pop++) {
                CurrentQ[pop] = Q[QPos (ind, pop)];
            }

            /*-------Do metropolis test of newq-------*/

            logdiff = 0.0;
            for (pop = 0; pop < MAXPOPS; pop++) {
                IndividualQ[pop] = TestQ[pop];
            }
            logdiff += Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                                Mapdistance, Phase,Phasemodel);

            for (pop = 0; pop < MAXPOPS; pop++) {
                IndividualQ[pop] = CurrentQ[pop];
            }
            logdiff -= Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                                Mapdistance, Phase,Phasemodel);

            randomnum = RandomReal (0.0, 1.0);
            if (randomnum < exp (logdiff)) {    /*accept */
                for (pop = 0; pop < MAXPOPS; pop++) {
                    Q[QPos (ind, pop)] = TestQ[pop];
                }
                numhits++;
            }
            /*for (pop=0;pop<MAXPOPS; pop++)
              printf("%1.3f %1.3f    ",CurrentQ[pop],TestQ[pop]);
              if (randomnum < exp(logdiff) ) printf("  Accepted ");
              printf("\n"); */

        }
    }
    if (REPORTHITRATE) {            /*does this once every UPDATEFREQ reps */
        if ((int) (rep - METROFREQ) / UPDATEFREQ < (int) rep / UPDATEFREQ) {
            printf ("Acceptance rate in UpdateQMetro %1.3f\n",
                    (float) numhits / NUMINDS);
        }
    }

    free (PriorQ);
    free (CurrentQ);
    free (TestQ);
    free (IndividualQ);
    free (RTransitProb);
}

/*---------------------------------------*/
/*void
  UpdateQRecombine (int *Geno, float *Q, int *Z, float *P,
  float *Alpha, int rep, struct IND *Individual,
  float *Mapdistance, float *R, float *Phase)
{
  int pop, ind;

  if (NOQS) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      Q[QPos (ind, pop)] = 1.0 / (float) MAXPOPS;
    }
  } else {
    UpdateQMetroRecombine (Geno, Q, Z, P, Alpha, rep,
                           Individual, Mapdistance, R, Phase,Phasemodel);
  }

}*/


/*-----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate LocPriors*/
void UpdateQCL (CLDict *clDict,int *Geno, int *PreGeno, float *Q, float *P, int *Z,
              float *Alpha, int rep,
              struct IND *Individual, float *UsePopProbs, int *Recessive, float *LocPrior,
              float * randomArr)
/*
 * update Q: proportion of ancest of each ind in each pop. Three
 * main options here.  One is update Q with admixture, one is update
 * Q with no admixture, and one is to use prior population info.
 * Even when USEPOPINFO is turned on, the other Q updates are still
 * called, in order to deal with any individuals who lack population
 * info.  The other functions are written to ignore individuals with
 * pop data when USEPOPINFO is on.
 */
{

    if (USEPOPINFO) {               /*update with prior population information */
        UpdateQWithPops (Geno, Q, P, Z, Alpha, rep, Individual, UsePopProbs,randomArr);
    }

    if (NOADMIX) {                  /*no admixture model */
        /* don't use ADMIBURNIN with LOCPRIOR models */
        if ((rep > ADMBURNIN) || (rep > BURNIN) || LOCPRIOR) {
            UpdateQNoAdmix (Geno, Q, P, Individual, LocPrior,randomArr);
        } else {
            UpdateQAdmixture (Q, Z, Alpha, Individual,randomArr);       /*initial burnin */
        }
    } else {
        /*admixture model */
        if (METROFREQ > 0 && rep%METROFREQ==0) {
            UpdateQMetroCL (clDict,Geno, PreGeno, Q, P, Alpha, rep, Individual, Recessive,
                          randomArr);
        } else {
            UpdateQAdmixtureCL (clDict,Q, Z, Alpha, Individual,randomArr);
        }
    }
}
