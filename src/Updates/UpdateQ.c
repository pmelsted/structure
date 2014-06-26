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

void UpdateQAdmixture (double *Q, int *Z, double *Alpha,
                       struct IND *Individual, double * randomArr);


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

/*----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate locprior*/
void UpdateQMetroCL (int *Geno, int *PreGeno, double *Q, double *P,
                   double *Alpha, int rep, struct IND *Individual, int *Recessive,
                   double * randomArr)

{
    int pop;
    double randomnum;
    int ind;
    int numhits = 0;

    double *TestQ;
    double *PriorQ1;
    double *logdiffs;

    RndDiscState randState[1];

    /* Removed:
     *   RECESSIVE ALLELES
     *   LOCPRIOR
     */


    TestQ = calloc (NUMINDS*MAXPOPS, sizeof (double));
    PriorQ1 = calloc (NUMINDS*MAXPOPS, sizeof (double));
    logdiffs = calloc(NUMINDS,sizeof(double));


    initRndDiscState(randState,randomArr,NUMINDS + NUMINDS*MAXRANDOM);
    rndDiscStateReset(randState, 0);
    /* ====== SAMPLE ======= */
    for(ind=0;ind < NUMINDS; ind++){
        for (pop = 0; pop < MAXPOPS; pop++) {
            PriorQ1[QPos(ind,pop)] = Alpha[pop];
        }
    }

    RDirichletDisc (PriorQ1, NUMINDS*MAXPOPS, TestQ, randState);

    /* ======== Calculate likelihood ====== */
    for (ind = 0; ind < NUMINDS; ind++) {
        logdiffs[ind] = 0.0;
        logdiffs[ind] += CalcLikeIndCL (Geno,TestQ,P, ind);
        logdiffs[ind] -= CalcLikeIndCL (Geno,Q, P, ind);
    }

    /* ========= Acceptance test ========= */
    for (ind = 0; ind < NUMINDS; ind++){
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            randomnum = rndDisc(randState);
            if (randomnum < exp (logdiffs[ind])) {    /*accept */
                for (pop = 0; pop < MAXPOPS; pop++) {
                    Q[QPos (ind, pop)] = TestQ[QPos(ind,pop)];
                }
                numhits++;
            }
        }
    }

    if (REPORTHITRATE) {           /*does this once every UPDATEFREQ reps */
        if ((int) (rep - METROFREQ) / UPDATEFREQ < (int) rep / UPDATEFREQ) {
            printf ("Acceptance rate in UpdateQMetro %1.3f\n",
                    (double) numhits / NUMINDS);
        }
    }

    free (PriorQ1);
    free (TestQ);
    free (logdiffs);
}



/*----------------------------------------*/
/*Melissa updated 7/12/07 to incorporate locprior*/
void UpdateQMetro (int *Geno, int *PreGeno, double *Q, double *P,
                   double *Alpha, int rep, struct IND *Individual, int *Recessive,
                   double * randomArr)
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
    double *PriorQ1;    /*[MAXPOPS]; */
    double *CurrentQ;             /*[MAXPOPS]; */
    double *TestQ;                /*[MAXPOPS]; */
    int pop;
    double logdiff;
    double randomnum;
    int ind;
    int numhits = 0;

    RndDiscState randState[1];
    /*  int i, ok; */

    /*  PriorQ1 = calloc (MAXPOPS, sizeof (double)); */
    CurrentQ = calloc (MAXPOPS, sizeof (double));
    TestQ = calloc (MAXPOPS, sizeof (double));

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
            RDirichletDisc (PriorQ1, MAXPOPS, TestQ, randState);

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
                    (double) numhits / NUMINDS);
        }
    }

    /*  free (PriorQ1); */
    free (CurrentQ);
    free (TestQ);
}


/*----------------------------------------*/
/*O(NUMINDS*MAXPOPS*LINES*NUMLOCI)*/
void
UpdateQNoAdmix (int *Geno, double *Q, double *P, struct IND *Individual,
                double *LocPrior,double * randomArr)
/*
 * Assign each individual to exactly one population according to the
 * conditional probabilities.
 */
{

    int ind, line, loc, pop;
    int allele;
    double *ProbsVector;          /*[MAXPOPS] */
    double sumlogs, sum;
    double runningtotal;
    double max=0.0, prob;
    int pickedpop;

    ProbsVector = calloc (MAXPOPS, sizeof (double));

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
void UpdateQAdmixture (double *Q, int *Z, double *Alpha,
                       struct IND *Individual,double * randomArr)
/*update Q: proportion of ancest of each ind in each pop. */
{
    int *NumLociPop;              /*[MAXPOPS] -- number of alleles from each pop (by ind) */
    double *Parameters;           /*[MAXPOPS] -- dirichlet parameters of posterior on Q */
    int ind, pop, loc;
    double *usealpha=NULL;

    NumLociPop = calloc (MAXPOPS, sizeof (int));
    Parameters = calloc (MAXPOPS, sizeof (double));
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
            RDirichlet (Parameters, MAXPOPS, Q + QPos (ind,
                        0));      /*generate dirichlet RVs */
        }
    }
    free (NumLociPop);
    free (Parameters);
}



/*-----------------------------------------*/
void
UpdateQWithPops (int *Geno, double *Q, double *P, int *Z, double *Alpha,
                 int rep, struct IND *Individual, double *UsePopProbs,
                 double * randomArr)
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
    double *Input, *Output;       /*both [MAXPOPS] */
    int ind, pop, loc, allele1, allele2, homepop;
    int pop1;
    int gen;
    double sumlogs;
    double sumsofar;
    double runningtotal;
    double *PostProbs;            /*[MAXPOPS][GENSBACK+1] */
    double *Priors;               /*[GENSBACK+1] */
    double p;
    double weights;
    double sumweights;
    double sumprobs;
    double rannum;
    double sqrtunder = sqrt (UNDERFLO);
    double like=0.0;              /*likelihood of current genotype */
    double maxprob;


    /*compute prior probs on each option */
    /*the correct-as-assigned option has weight MIGRPRIOR; the remaining
      probability is split among the other options */

    Input = calloc (MAXPOPS, sizeof (double));
    Output = calloc (MAXPOPS, sizeof (double));
    PostProbs = calloc (MAXPOPS * (GENSBACK + 1), sizeof (double));
    Priors = calloc (GENSBACK + 1, sizeof (double));
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
                                   + ((double) 1.0 - p) * P[PPos (loc, homepop, allele1)];
                        } else if (allele2 != MISSING) {
                            like = p * P[PPos (loc, homepop, allele2)]
                                   + ((double) 1.0 - p) * P[PPos (loc, homepop, allele2)];
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
void UpdateQ (int *Geno, int *PreGeno, double *Q, double *P, int *Z,
              double *Alpha, int rep,
              struct IND *Individual, double *UsePopProbs, int *Recessive, double *LocPrior,
              double * randomArr)
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
            UpdateQMetroCL (Geno, PreGeno, Q, P, Alpha, rep, Individual, Recessive,
                          randomArr);
        } else {
            UpdateQAdmixture (Q, Z, Alpha, Individual,randomArr);
        }
    }
}

void
UpdateQMetroRecombine (int *Geno, double *Q, int *Z, double *P,
                       double *Alpha, int rep, struct IND *Individual,
                       double *Mapdistance, double *R,
                       double *Phase,int *Phasemodel,
                       double * randomArr)
/*
 * This function does the same job as UpdateQMetro in the case
 * when there is recombination.
 * the code is very similar, but a new routine is nevertheless a good idea
 */
{
    double *PriorQ;               /*[MAXPOPS]; */
    double *CurrentQ;             /*[MAXPOPS]; */
    double *TestQ;                /*[MAXPOPS]; */
    int pop;
    double logdiff;
    double randomnum;
    int ind;
    int numhits = 0;
    /*  int i, ok; */
    double *RTransitProb;
    double *IndividualQ;

    if (PHASED) {
        RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
    } else {
        RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
    }

    PriorQ = calloc (MAXPOPS, sizeof (double));
    CurrentQ = calloc (MAXPOPS, sizeof (double));
    TestQ = calloc (MAXPOPS, sizeof (double));
    IndividualQ = calloc (MAXPOPS, sizeof (double));

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
                    Q[QPos (ind, pop)] = (double) 1 / MAXPOPS;
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
                    (double) numhits / NUMINDS);
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
  UpdateQRecombine (int *Geno, double *Q, int *Z, double *P,
  double *Alpha, int rep, struct IND *Individual,
  double *Mapdistance, double *R, double *Phase)
{
  int pop, ind;

  if (NOQS) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      Q[QPos (ind, pop)] = 1.0 / (double) MAXPOPS;
    }
  } else {
    UpdateQMetroRecombine (Geno, Q, Z, P, Alpha, rep,
                           Individual, Mapdistance, R, Phase,Phasemodel);
  }

}*/

