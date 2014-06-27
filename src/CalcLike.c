#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "Kernels.h"

/*-----------------------------------------*/
double
CalcLikeIndRecessive(int *Geno, int *PreGeno, double *AncVector, double *P,
                     int ind, int *Recessive)
{
    /*returns log(likelihood) of Data for one individual:  log[ P(Data|p,q) ].
      This version is used only for diploid individuals when there is genotypic
      uncertainty (recessive model or inbreeding [in future])

      note code overlap with CalcLikeRecessive; make any updates to both functions

      Notice use of AncVector (Q for current individual only), not full Q)
    */

    double runningtotal = 1;
    double loglike = 0;
    double term, sum1, sum2;
    int allele1, allele2;
    int loc, pop;
    /*  int line; */

    for (loc = 0; loc < NUMLOCI; loc++) {
        allele1 = PreGeno[GenPos (ind, 0, loc)];
        allele2 = PreGeno[GenPos (ind, 1, loc)];

        if (allele1 != MISSING && allele2 != MISSING) { /* no missing data */
            sum1 = 0.0;
            sum2 = 0.0;
            for (pop = 0; pop < MAXPOPS; pop++) {
                /* summing over possible genotypes and possible Z */
                sum1 += AncVector[pop] * P[PPos (loc, pop, allele1)];
                if (Recessive[loc] != MISSING && Recessive[loc] != allele1
                        && allele1 == allele2) { /* bug fixed 05072007 */
                    sum2 += AncVector[pop] * (P[PPos (loc, pop, allele1)] + P[PPos(loc, pop,
                                              Recessive[loc])]);
                } else {
                    sum2 += AncVector[pop] * P[PPos (loc, pop, allele2)];
                }
            }
            term = sum1 * sum2;
        } else if (allele1!=MISSING)  { /* one allele missing */
            term=0.0;
            for (pop=0; pop<MAXPOPS; pop++) {
                term += AncVector[pop] * P[PPos (loc, pop, allele1)];
            }
        } else if (allele2!=MISSING) {
            term=0.0;
            for (pop=0; pop<MAXPOPS; pop++) {
                term += AncVector[pop] * P[PPos (loc, pop, allele2)];
            }
        } else {
            term=1.0;  /* no data */
        }

        runningtotal *= term;

        if (runningtotal <
                UNDERFLO) {        /*this is to avoid having to take logs all the time */
            loglike += log (runningtotal);
            runningtotal = 1;
        }
    }

    loglike += log (runningtotal);
    return loglike;
}

/*-----------------------------------------*/
double CalcLikeInd (int *Geno, int *PreGeno, double *AncVector, double *P,
                    int ind, int *Recessive)
{
    /*returns log(likelihood) of Data for one individual:  log[ P(Data|p,q) ]
      See notes 19 March 99 */

    /*This likelihood is used for the metropolis update of Q*/

    /* when there is genotypic ambiguity (eg in the recessive model) the likelihood
       is computed in one of two ways: (1) for diploids it sums over possible
       genotypes, and (2) for other ploidy it is computed based on the current
       imputed genotypes

       note code overlap with CalcLike, CalcLikeIndRecessive;
       make any updates to all functions

       Notice use of AncVector (Q for current individual only), not full Q)
    */

    double runningtotal = 1;
    double loglike = 0;
    double term;
    int allele;
    int line, loc, pop;
    double sqrtunder = sqrt (UNDERFLO);

    if (LINES==2 && RECESSIVEALLELES) {
        loglike = CalcLikeIndRecessive(Geno, PreGeno, AncVector, P, ind, Recessive);
    } else {
        for (line = 0; line < LINES; line++) {
            for (loc = 0; loc < NUMLOCI; loc++) {
                allele = Geno[GenPos (ind, line, loc)];
                if (allele != MISSING) {
                    term = 0.0;
                    for (pop = 0; pop < MAXPOPS; pop++) {
                        term += AncVector[pop] * P[PPos (loc, pop, allele)];
                    }

                    if (term > sqrtunder) {
                        runningtotal *= term;
                    } else {
                        runningtotal *= sqrtunder;
                    }

                    if (runningtotal <
                            sqrtunder) { /*this is to avoid having to take logs all the time */
                        if (runningtotal == 0.0) {
                            printf ("Error in CalcLikeInd\n");
                        }
                        loglike += log (runningtotal);
                        runningtotal = 1;
                    }
                }
            }
        }
        loglike += log (runningtotal);
    }
    return loglike;
}


double CalcLikeIndDiffCL (int *Geno,  double *TestQ, double *Q, double *P, int ind)
{

    double logdiff;
    double termPlus, termMinus,logPlus = 0.0,logMinus = 0.0;
    int allele;
    int line, loc, pop;

    for (line = 0; line < LINES; line++) {
        for (loc = 0; loc < NUMLOCI; loc++) {
            allele = Geno[GenPos (ind, line, loc)];
            if (allele != MISSING) {
                termPlus = 0.0;
                termMinus = 0.0;
                for (pop = 0; pop < MAXPOPS; pop++) {
                    termPlus += TestQ[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                    termMinus += Q[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                }
                logPlus += log(termPlus);
                logMinus += log(termMinus);
            }
        }
    }
    logdiff = logPlus - logMinus;
    return logdiff;
}

double CalcLikeIndCL (int *Geno, double *Q, double *P, int ind)
{

    double term;
    double logterm = 0.0;
    int allele;
    int line, loc, pop;

    for (loc = 0; loc < NUMLOCI; loc++) {
        for (line = 0; line < LINES; line++) {
            allele = Geno[GenPos (ind, line, loc)];
            if (allele != MISSING) {
                term = 0.0;
                for (pop = 0; pop < MAXPOPS; pop++) {
                    term += Q[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                }
                logterm += log(term);
            }
        }
    }
    return logterm;
}

void reduceLogdiffsCL(double *logterms, double *logdiffs){
    double logterm;
    int ind, loc;
    for (ind =0; ind < NUMINDS; ind++){
        logterm = 0.0;
        for (loc = 0; loc < NUMLOCI; loc++) {
            logterm += logterms[ind*NUMLOCI + loc];
        }
        logdiffs[ind] = logterm;
    }
}

void CalcLogdiffsCL(CLDict *clDict,int *Geno,double *TestQ, double *Q, double *P, double *logdiffs)
{

    double termP,termM;
    double logterm;
    int allele;
    int line, loc, pop,ind;
    double * logdiffsnoncl;
    double *logterms;
    size_t global[2];
    global[0] = NUMINDS;
    global[1] = NUMLOCI;

    logterms = calloc(NUMINDS*NUMLOCI,sizeof(double));

    if(DEBUGCOMPARE){
        logdiffsnoncl = calloc(NUMINDS,sizeof(double));
        for (ind =0; ind < NUMINDS; ind++){
            logterm = 0.0;
            for (loc = 0; loc < NUMLOCI; loc++) {
                logterms[ind*NUMLOCI + loc] = 0.0;
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
            }
            logdiffsnoncl[ind] = logterm;
        }
    }


    writeBuffer(clDict,Q,sizeof(double) * QSIZE,QCL,"Q");
    writeBuffer(clDict,TestQ,sizeof(double) * QSIZE,TESTQCL,"TestQ");
    /* Already up to date on GPU */
    /*
    writeBuffer(clDict,P,sizeof(double) * PSIZE,PCL,"P");
    writeBuffer(clDict,Geno,sizeof(int) * GENOSIZE,GENOCL,"GENO");
    */


    runKernel(clDict,mapLogDiffsKernel,2,global,"mapLogDiffs");
    readBuffer(clDict,logterms,sizeof(double) * NUMINDS*NUMLOCI,LOGTERMSCL,"Logterms");
    reduceLogdiffsCL(logterms,logdiffs);

    if (DEBUGCOMPARE){
        for(ind = 0; ind < NUMINDS; ind++){
            if(fabs(logdiffs[ind]-logdiffsnoncl[ind]) > 10e-6){
                    printf("C %f G %f\n ",
                           logdiffs[ind],
                           logdiffsnoncl[ind]
                        );
                    ReleaseCLDict(clDict);
                    exit(EXIT_FAILURE);
            }
        }
    }



    /*for(ind = 0; ind < NUMINDS; ind++){
        for(loc=0;loc < NUMLOCI; loc++){
            if(fabs(logtermscl[ind]-logterms[ind]) > 10e-8){
                printf("%f, %f, %d,%d\n", logtermscl[ind*NUMLOCI + loc],logterms[ind*NUMLOCI + loc], ind,loc);
            }
        }
    }*/

    free(logterms);

    if(DEBUGCOMPARE){
        free(logdiffsnoncl);
    }

}



double CalcLikeIndDiff (int *Geno, int *PreGeno, double *AncVectorPlus,
                        double *AncVectorMinus, double *P, int ind, int *Recessive)
{
    /*returns log(likelihood) of Data for one individual:  log[ P(Data|p,q) ]
      See notes 19 March 99 */

    /*This likelihood is used for the metropolis update of Q*/

    /* when there is genotypic ambiguity (eg in the recessive model) the likelihood
       is computed in one of two ways: (1) for diploids it sums over possible
       genotypes, and (2) for other ploidy it is computed based on the current
       imputed genotypes

       note code overlap with CalcLike, CalcLikeIndRecessive;
       make any updates to all functions

       Notice use of AncVector (Q for current individual only), not full Q)
    */

    /*double runningtotalPlus = 1, runningtotalMinus =1;*/
    double logdiff = 0;
    double termPlus, termMinus,logPlus = 0.0,logMinus = 0.0;
    int allele;
    int line, loc, pop;
    /*double sqrtunder = sqrt (UNDERFLO);*/

    for (line = 0; line < LINES; line++) {
        for (loc = 0; loc < NUMLOCI; loc++) {
            allele = Geno[GenPos (ind, line, loc)];
            if (allele != MISSING) {
                termPlus = 0.0;
                termMinus = 0.0;
                for (pop = 0; pop < MAXPOPS; pop++) {
                    termPlus += AncVectorPlus[pop] * P[PPos (loc, pop, allele)];
                    termMinus += AncVectorMinus[pop] * P[PPos (loc, pop, allele)];
                }
                logPlus += log(termPlus);
                logMinus += log(termMinus);

                /*runningtotalPlus *= termPlus > sqrtunder ? termPlus : sqrtunder;
                runningtotalMinus *= termMinus > sqrtunder ? termMinus : sqrtunder;

                if (runningtotalPlus <
                        sqrtunder) { [>this is to avoid having to take logs all the time <]
                    if (runningtotalPlus == 0.0) {
                        printf ("Error in CalcLikeInd\n");
                    }
                    logdiff += log (runningtotalPlus);
                    runningtotalPlus = 1;
                }

                if (runningtotalMinus <
                        sqrtunder) { [>this is to avoid having to take logs all the time <]
                    if (runningtotalMinus == 0.0) {
                        printf ("Error in CalcLikeInd\n");
                    }
                    logdiff -= log (runningtotalMinus);
                    runningtotalMinus = 1;
                }*/
            }
        }
    }
    logdiff = logPlus - logMinus;
    return logdiff;
}
/*-----------------------------------------*/

double CalcLikeRecessive(int *Geno, int *PreGeno, double *Q, double *P,
                         int *Recessive, double *sumindlike, double *indlike_norm)
/* calculate log likelihood for recessive model.  Only coded for diploids!
 *
 * note code overlap with CalcLikeInd, CalcLikeIndRecessive, CalcLike;
 * make any updates to all functions */

/*7/12/07: Melissa modified so that it calls CalcLikeIndRecessive
  and there is no more code overlap.  This is more efficient and better
  coding practice, and also helps with calculation of DIC (which
  requires storage of individual likelihoods)*/

{


    double loglike = 0;
    int ind, pop;
    double *AncVector = malloc(MAXPOPS*sizeof(double));
    double indlike;
    if (AncVector==NULL) {
        printf("Error allocating Ancvetor in CalcLikeRecessive: out of memory?\n");
        exit(-1);
    }

    for (ind = 0; ind < NUMINDS; ind++) {
        for (pop=0; pop<MAXPOPS; pop++) {
            AncVector[pop] = Q[QPos(ind, pop)];
        }
        indlike = CalcLikeIndRecessive(Geno, PreGeno, AncVector, P, ind,
                                       Recessive);
        if (sumindlike!=NULL) {
            sumindlike[ind] += exp(indlike-indlike_norm[ind]);
        }
        loglike += indlike;
    }
    free(AncVector);
    return loglike;
}
/*-----------------------------------------*/
double CalcLike (int *Geno, int *PreGeno, double *Q, double *P, int *Recessive,
                 double *sumindlike, double *indlike_norm)
{
    /*returns log(likelihood) of Data:  log[ P(Data|p,q) ]
      See notes 19 March 99 */

    /* note code overlap with CalcLikeInd, CalcLikeIndRecessive, CalcLikeRecessive;
     * make any updates to all functions */
    /*Melissa modified 7/12/07 so it calls CalcLikeInd rather than repeats
      the same code.  This helps with DIC calculation...*/

    double loglike = 0, indlike, *AncVector;
    int ind, pop;

    if (LINES == 2 && RECESSIVEALLELES) {
        loglike = CalcLikeRecessive(Geno, PreGeno, Q, P, Recessive, sumindlike,
                                    indlike_norm);
    } else {
        AncVector = malloc(MAXPOPS*sizeof(double));
        if (AncVector==NULL) {
            printf("Error allocating Ancvetor in CalcLikeRecessive: out of memory?\n");
            exit(-1);
        }
        for (ind = 0; ind < NUMINDS; ind++) {
            for (pop=0; pop<MAXPOPS; pop++) {
                AncVector[pop] = Q[QPos(ind, pop)];
            }
            indlike = CalcLikeInd(Geno, PreGeno, AncVector, P, ind, Recessive);
            if (sumindlike!=NULL) {
                if (indlike_norm[ind]==0.0) {
                    indlike_norm[ind] = indlike;
                }
                sumindlike[ind] += exp(indlike-indlike_norm[ind]);
            }
            loglike += indlike;
        }
        free(AncVector);
    }
    return loglike;
}
