    /*Part of structure.c.

  This bit of the program is involved in collecting data (DataCollection)
  and printing results (OutPutResults). */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "structure.h"
#include "params.h"
#include "mymath.h"
#include "CalcLike.h"
#include "Kernels.h"

void UpdateSums (float *Q, float *QSum, int *Z, float *P, float *PSum,
                 float *Fst, float *FstSum, int *NumAlleles,
                 int *AncestDist, float *Epsilon, float *SumEpsilon,
                 float *lambda, float *sumlambda,
                 float *LocPrior, float *sumLocPrior, int LocPriorLen);
float CalcLike (int *Geno, int *PreGeno, float *Q, float *P, int *Recessive,
                 float *sumindlike, float *indlike_norm);
float EstLogProb (float sumlikes, float sumsqlikes, int reps);
float KLDiv (int pop1, int pop2, float *P, float *LogP, int *NumAlleles,
              int reps);
void PrintKLD (FILE * file, float *P, float *LogP, int *NumAlleles, int reps,
               int format);
void PrintNET (FILE * file, float *P, int *NumAlleles, int reps, int format);

void PrintBanner (int rep, float *Alpha, float *Fst, float like,
                  float *lambda);
void PrintMainParams (FILE * file, int rep, int argc, char *argv[]);
void PrintQ (FILE * file, int *Geno, int rep, float *QSum,
             struct IND *Individual,
             int *AncestDist, float *UsePopProbs,float *sumR);
void PrintP (FILE * file, int rep, int *Geno, float *PSum, int *Translation,
             int *NumAlleles, float *SumEpsilon, char *Markername);
void PrintMembership (FILE * file, float *QSum, struct IND *Individual);
void PrintSequences (FILE * file, int *Geno, char *Markername, int rep,
                     int *Translation);
void PrintSums (FILE * file, int rep, float sumlikes,
                float sumsqlikes, float *FstSum, float *sumAlpha, float *sumlambda,
                float *sumR, float *varR, struct IND *Individual,
                float *sumLocPriors, int LocPriorLen, float DIC);
void PrintGeneName(FILE * file, int loc, char *Markername);
int EqualGeneNames(int loc1,int loc2,char *Markername);


/*=================================================*/
void
DataCollection (int *Geno, int *PreGeno,
                float *Q, float *QSum, int *Z, int *Z1,
                float *P, float *PSum,
                float *Fst, float *FstSum, int *NumAlleles,
                int *AncestDist, float *Alpha, float *sumAlpha,
                float *sumR, float *varR, float *like,
                float *sumlikes, float *sumsqlikes, float *R,
                float *Epsilon, float *SumEpsilon, float recomblikelihood,
                float *lambda, float *sumlambda, int *Recessive,
                float *LocPrior, float *sumLocPrior, int LocPriorLen,
                float *sumindlikes, float *indlikes_norm, int rep)
{
    int ind, pop, loc, pos;
    UpdateSums (Q, QSum, Z, P, PSum, Fst, FstSum, NumAlleles, AncestDist,
                Epsilon, SumEpsilon, lambda, sumlambda,
                LocPrior, sumLocPrior, LocPriorLen);
    if (LINKAGE) {
        for (ind = 0; ind < NUMINDS; ind++) {
            sumR[ind] += R[ind];
            varR[ind] += R[ind] * R[ind];
        }
    }

    if (LOCPRIOR && NOADMIX==0) {
        for (pop=0; pop<MAXPOPS; pop++)
            for (loc=0; loc<=NUMLOCATIONS; loc++) {
                pos = AlphaPos(loc, pop);
                sumAlpha[pos] += Alpha[pos];
            }
    } else if (!(NOADMIX) && (!(NOALPHA))) {
        for (pop = 0; pop < MAXPOPS; pop++) {
            sumAlpha[pop] += Alpha[pop];
        }
    }


    if (COMPUTEPROB) {
        if (LINKAGE) {
            *like = recomblikelihood;
        }

        if (rep < BURNIN) {
            *like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
        } else {
            *like = CalcLike (Geno, PreGeno, Q, P, Recessive,
                              sumindlikes, indlikes_norm);
        }
        *sumlikes += *like;
        *sumsqlikes += (*like) * (*like);
    }
    /*printf("%f %f %f\n", *like, *sumlikes, *sumsqlikes); */
}

void
DataCollectionCL (CLDict *clDict,int *Geno, int *PreGeno,
                float *Q, float *QSum, int *Z, int *Z1,
                float *P, float *PSum,
                float *Fst, float *FstSum, int *NumAlleles,
                int *AncestDist, float *Alpha, float *sumAlpha,
                float *sumR, float *varR, float *like,
                float *sumlikes, float *sumsqlikes, float *R,
                float *Epsilon, float *SumEpsilon, float recomblikelihood,
                float *lambda, float *sumlambda, int *Recessive,
                float *LocPrior, float *sumLocPrior, int LocPriorLen,
                float *sumindlikes, float *indlikes_norm, int rep)
{
    int ind, pop, loc, pos;
    int i;
    int usesumindlikes;
    float gpulike[1];
    size_t global[2];

    if (LOCPRIOR){
        for (i=0; i<LocPriorLen; i++) {
            sumLocPrior[i] += LocPrior[i];
        }
    }

    if (LINKAGE) {
        for (ind = 0; ind < NUMINDS; ind++) {
            sumR[ind] += R[ind];
            varR[ind] += R[ind] * R[ind];
        }
    }

    global[0] = MAXPOPS;
    runKernel(clDict,DataCollectPopKernel,1,global,"DataCollectionPop");
    global[1] = NUMINDS;
    runKernel(clDict,DataCollectIndKernel,2,global,"DataCollectionInd");
    global[1] = NUMLOCI;
    runKernel(clDict,DataCollectLocKernel,2,global,"DataCollectionLoc");

    if (COMPUTEPROB) {
        /*if (LINKAGE) {*/
            /**like = recomblikelihood;*/
        /*}*/

        if (rep == BURNIN) {
            usesumindlikes = 1;
            setKernelArgExplicit(clDict,CalcLikeKernel,sizeof(int),&usesumindlikes,3);
        }
        global[1] = NUMINDS;
        /*global[0] = (32 < NUMLOCI) ? 32 : NUMLOCI;*/
        /*global[0] = ( 128 < NUMLOCI) ? 128 : NUMLOCI;*/
        global[0] = pow(2,(int) (log(NUMLOCI)/log(2)));
        runKernel(clDict,MapReduceLogLikeKernel,2,global,"LogLike");

        global[0] = NUMINDS;
        runKernel(clDict,CalcLikeKernel,1,global,"CalcLike");
        global[0] = 1;
        runKernel(clDict,ComputeProbFinishKernel,1,global,"ComputeProbFinish");
        global[1] = NUMLOCI;

        /*readBuffer(clDict,gpulike, sizeof(float),LIKECL, "like");*/

        /*if (rep < BURNIN) {*/
            /**like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);*/
        /*} else {*/

            /**like = CalcLike (Geno, PreGeno, Q, P, Recessive,*/
                              /*sumindlikes, indlikes_norm);*/
        /*}*/
        /*if (rep % 100 == 0){*/
            /*printf("like %f %f ",*like,*gpulike);*/
            /*if(fabs(*like - *gpulike) > 10e-10){*/
                /*printf("DIFF!");*/
            /*}*/
            /*printf("\n");*/
        /*}*/
        /**sumlikes += *like;*/
        /**sumsqlikes += (*like) * (*like);*/
    }
    /*printf("%f %f %f\n", *like, *sumlikes, *sumsqlikes); */
}

/*---------------------------------------------------*/
void
PrintLike (float like, int rep, int *Geno, int *PreGeno, float *Q,
           float *P,float recomblikelihood,
           int *Recessive)
{
    if (rep + 1 > BURNIN) {        /*already calculated */
        printf ("%6.0f\n", like);
    } else {
        if (LINKAGE) {
            printf("%6.0f\n",recomblikelihood);
        } else {
            printf ("%6.0f\n", CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL));
        }
    }
}


/*---------------------------------------------------*/
float
CalculateRAverage (float *R)
{
    int ind;
    float temp;
    temp = 0.0;
    for (ind = 0; ind < NUMINDS; ind++) {
        temp += R[ind];
    }
    return temp / (float) NUMINDS;
}

/*----------------------------------------------------*/
void
PrintUpdate (int rep, int *Geno, int *PreGeno, float *Alpha, float *Fst,
             float *P, float *Q,
             float like, float sumlikes, float sumsqlikes, int *NumAlleles,
             float *R, float *lambda, struct IND *Individual,
             float recomblikelihood, int *Recessive,
             float *LocPrior, int LocPriorLen)
{
    /*print a bunch of stuff to screen during run: rep, alpha, f, KLD, likelihood...
      Also occasionally print a header banner to define the variables. */

    float logprob=0.0;
    /*  int i;
     *  int printalign; */
    int pop;

    if ((rep < BURNIN + UPDATEFREQ) && (rep > BURNIN)) {
        printf ("\nBURNIN completed");
    }

    if (((NOADMIX) && (ADMBURNIN > 0)) && !LOCPRIOR ) {
        if ((rep < ADMBURNIN + UPDATEFREQ) && (rep >= ADMBURNIN)) {
            printf ("\nAdmixture Burnin complete.  Current alpha");
            if (POPALPHAS) {
                printf ("s = ");
                for (pop = 0; pop < MAXPOPS; pop++) {
                    printf ("%1.3f ", Alpha[pop]);
                }
            } else {
                printf (" = %1.3f ", Alpha[0]);
            }
            printf ("\n\n");
        }
    }

    if ((LINKAGE) && (ADMBURNIN> 0)) {
        if ((rep < ADMBURNIN+UPDATEFREQ) && (rep >= ADMBURNIN)) {
            printf ("\nNo recombination Burnin complete.  Current rec = %1.3f",
                    CalculateRAverage (R));
            PrintBanner (rep, Alpha, Fst, like, lambda);
        }
    }

    /*calculate some stuff */

    if (LINKAGE) {
        if (rep <= ADMBURNIN) {
            like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
        } else {
            like = recomblikelihood;
        }

        if (rep >= BURNIN + 2) { /* +2 because need 2 observations for variance*/
            logprob = EstLogProb (sumlikes, sumsqlikes, rep - BURNIN);
        } else if (COMPUTEPROB) { /*not linkage model*/
            if (rep <= BURNIN + 2) {
                like = CalcLike (Geno, PreGeno, Q, P, Recessive, NULL, NULL);
            } else {
                logprob = EstLogProb (sumlikes, sumsqlikes, rep - BURNIN);
            }
        }
    }

    /*possibly print banner to indicate what the numbers are */
    if ((rep == UPDATEFREQ) ||
            ((rep < BURNIN) && ((rep % (BANNERFREQ * UPDATEFREQ)) == 0)) ||
            ((rep > BURNIN) && (((rep - BURNIN) % (BANNERFREQ * UPDATEFREQ)) == 0))) {
        if (rep != NUMREPS + BURNIN) {        /*don't bother for last line of output */
            if ((rep > 1) && (PRINTQSUM) && (MAXPOPS > 1)) {
                PrintMembership (stdout, Q, Individual);
            }
            PrintBanner (rep, Alpha, Fst, like, lambda);
        }
    }

    /*print current values to screen */
    printf ("%5d:    ", rep);

    if (LINKAGE && rep >= ADMBURNIN) {
        if (!INDIVIDUALR) {
            printf ("%1.09f  ", R[0]);
        } else {
            printf ("%1.09f   ", CalculateRAverage (R));
        }
    }

    if (PRINTLAMBDA) {
        if (POPSPECIFICLAMBDA) {
            for (pop=0; pop<MAXPOPS; pop++) {
                printf ("%1.2f    ", lambda[pop]);
            }
        } else {
            printf ("%1.2f    ", lambda[0]);
        }
    }

    if ((!(NOADMIX)) && (!(NOALPHA))) {
        if (POPALPHAS) {
            for (pop = 0; pop < MAXPOPS; pop++) {
                printf ("%1.3f  ", Alpha[pop]);
                if (pop > 8) {
                    printf (" "); /*extra space for number */
                }
            }
        } else {
            printf ("%1.3f  ", Alpha[0]);
        }
    }

    if (FREQSCORR) {
        printf ("  ");
        if (ONEFST) {
            printf ("%1.3f ", Fst[0]);
        } else {
            for (pop = 0; pop < MAXPOPS; pop++) {
                printf ("%1.3f ", Fst[pop]);
            }
        }
        printf ("   ");
    } else {
        printf ("  ");
    }

    if (LOCPRIOR) {
        printf ("%1.3f ", LocPrior[0]);
        printf ("   ");
    }

    /*currently it only net distances not KLD ones */
    if (PRINTKLD || PRINTNET) {
        PrintNET (stdout, P, NumAlleles, 1, 0);
    }

    if (COMPUTEPROB) {


        /*put correct # of spaces in */
        /*if (((int) log10(like)) < 1) printalign = 4;
          else printalign = 4 + ((int) log10(like));
          for (i=printalign; i<8; i++)
          printf(" "); */
        if (rep > BURNIN + 2) {
            printf ("  %.0f  ", like);
            printf ("  %.0f ", logprob);
        } else {
            printf ("  --  ");
        }
    }

    printf ("\n");

    if (rep == BURNIN) {
        printf ("\nBURNIN completed");
        PrintBanner (rep, Alpha, Fst, like, lambda);
    }
    fflush(stdout);

}
/*----------------------------------------------------*/
void
PrintBanner (int rep, float *Alpha, float *Fst, float like, float *lambda)
/*print banner to screen during run giving variable names */
{
    int i, j, k;
    int pop;

    printf ("\n");
    for (i = 4; i < ((int) log10 (rep)); i++) {
        printf (" ");
    }
    printf (" Rep#:   ");

    if ((LINKAGE) && (rep >= ADMBURNIN)) {
        printf (" r           ");
    }

    if (PRINTLAMBDA) {
        if (POPSPECIFICLAMBDA) {
            for (pop=0; pop<MAXPOPS; pop++) {
                printf("Lambda%d ",pop+1);
            }
        } else {
            printf ("Lambda  ");
        }
    }

    if (((!(NOADMIX)) && (!(NOALPHA)))) {
        printf (" ");
        if (POPALPHAS) {
            for (pop = 0; pop < MAXPOPS; pop++) {
                printf ("Alpha%d ", pop + 1);
            }
        } else {
            printf ("Alpha  ");
        }
        /*for (i = 7; i < ((int) log10 (Alpha[0])); i++) {
          printf (" ");
        }*/
    }

    if (FREQSCORR) {
        printf ("   ");
        if (ONEFST) {
            printf ("Fst   ");
        } else
            for (pop = 0; pop < MAXPOPS; pop++) {
                printf (" F%d   ", pop + 1);
            }
        printf (" ");
    } else {
        printf (" ");
    }

    if (LOCPRIOR) {
        printf ("  r     ");
    }

    if (PRINTKLD || PRINTNET) {
        for (j = 0; j < MAXPOPS - 1; j++)
            for (k = j + 1; k < MAXPOPS; k++) {
                printf (" D%d,%d ", j + 1, k + 1);
            }
    }

    if (COMPUTEPROB) {
        printf ("   Ln Like ");
        for (i = 8; i < ((int) log10 (rep)); i++) {
            printf (" ");
        }

        if (rep >= BURNIN) {
            printf (" Est Ln P(D)");
        }
    }
    printf ("\n");


}

/*----------------------------------------------------*/
float
EstLogProb (float sumlikes, float sumsqlikes, int reps)
/*returns the current estimated Prob of Data.  Reps is the
  number of reps, not including burnin */
{
    float mean = sumlikes / reps;
    float var = SampleVar (sumsqlikes, sumlikes, reps);

    return (mean - var / 2.0);

}

/*-------------------------------------------------------*/
float
NETiv (int pop1, int pop2, float *P, int *NumAlleles, int reps)
{
    /* This function returns the current estimated average net nucleotide
       distance between the allele frequencies in populations 1 and 2.
       Here reps is the number of reps over which the P is an average, rather
       than the number of reps since the start of the program, as elsewhere */
    float sum, d1, d2,d12, norm;
    int loc,allele;
    norm = (float) reps;
    norm *= norm;
    sum=0.0;
    for (loc=0; loc<NUMLOCI; loc++) {
        d1=0.0;
        d2=0.0;
        d12=0.0;
        for (allele=0; allele<NumAlleles[loc]; allele++) {
            d1+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop1,allele)]/norm;
            d2+=P[PPos(loc,pop2,allele)]*P[PPos(loc,pop2,allele)]/norm;
            d12+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop2,allele)]/norm;
        }
        sum+= 0.5*(d1+d2)-d12;
    }
    return sum/NUMLOCI;
}


/*-------------------------------------------------------*/
float
GROSSiv (int pop1, int pop2, float *P, int *NumAlleles, int reps)
{
    /* This function returns the current estimated average gross nucleotide
       distance between the allele frequencies in populations 1 and 2.
       Here reps is the number of reps over which the P is an average, rather
       than the number of reps since the start of the program, as elsewhere */
    float sum, d12, norm;
    int loc,allele;
    norm = (float) reps;
    norm *= norm;
    sum=0.0;
    for (loc=0; loc<NUMLOCI; loc++) {
        d12=0.0;
        for (allele=0; allele<NumAlleles[loc]; allele++) {
            d12+=P[PPos(loc,pop1,allele)]*P[PPos(loc,pop2,allele)]/norm;
        }
        sum+= 1.0-d12;
    }
    return sum/NUMLOCI;
}

/*----------------------------------------------------*/
float
KLDiv (int pop1, int pop2, float *P, float *LogP, int *NumAlleles, int reps)
/*This function returns the current (average) estimated
  Kullback-Leibler divergence between the allele frequencies in
  pops 1 and 2.  Here reps is the number of reps over which the P
  is an average, rather than the number of reps since the start of
  the program, as elsewhere. */
{
    float sum = 0.0;
    int allele, loc;

    for (loc = 0; loc < NUMLOCI; loc++)
        for (allele = 0; allele < NumAlleles[loc]; allele++)
            sum += ((float) P[PPos (loc, pop1, allele)] / reps)
                   * log (P[PPos (loc, pop1, allele)] / P[PPos (loc, pop2, allele)]);

    return sum / NUMLOCI;

}

/*----------------------------------------------------*/
void
PrintNET (FILE * file, float *P, int *NumAlleles, int reps, int format)
/*This function prints the current (average) estimated
  Net-nucleotide divergence between the allele frequencies in the
  different populations, in two formats.  Format 0 prints these in a
  row, averaging D_ij with D_ji, while Format 1 prints the full table */
{
    int i, j;
    /* int k;
     *  float div; */

    if (format == 0) {             /*print in line */
        for (i = 0; i < MAXPOPS - 1; i++) {
            for (j = i + 1; j < MAXPOPS; j++) {
                fprintf (file, "%1.3f ", NETiv(i,j,P,NumAlleles,reps));
            }
        }
    } else {
        /*print 2-D table */
        fprintf (file,
                 "\nAllele-freq. divergence among pops (Net nucleotide distance)");
        if (reps > 1) {
            fprintf (file, ",\ncomputed using point estimates of P");
        }
        fprintf (file, ".\n\n");

        fprintf (file, "     ");
        for (j = 0; j < MAXPOPS; j++) {
            fprintf (file, "%2d      ", j + 1);
        }
        fprintf (file, "\n");

        for (i = 0; i < MAXPOPS; i++) {
            fprintf (file, "%2d   ", i + 1);
            for (j = 0; j < MAXPOPS; j++) {
                if (i == j) {
                    fprintf (file, "   -    ");
                } else {
                    fprintf (file, "%1.4f  ", NETiv (i, j, P, NumAlleles, reps));
                }
            }
            fprintf (file, "\n");
        }

        /* word change population -> cluster, William 03/27/07 */
        fprintf(file,
                "\nAverage distances (expected heterozygosity) between individuals in same cluster:\n");
        for (i=0; i<MAXPOPS; i++) {
            fprintf(file,"cluster %2d  : %1.4f \n",i+1,GROSSiv(i,i,P,NumAlleles,reps));
        }
        fprintf(file,"\n");
    }
}


/*----------------------------------------------------*/
void
PrintKLD (FILE * file, float *P, float *LogP, int *NumAlleles, int reps,
          int format)
/*This function prints the current (average) estimated
  Kullback-Leibler divergence between the allele frequencies in the
  different populations, in two formats.  Format 0 prints these in a
  row, averaging D_ij with D_ji, while Format 1 prints the full table */
{
    int i, j;
    int k;
    float div;

    if (format == 0) {            /*print in line */
        for (i = 0; i < MAXPOPS - 1; i++)
            for (j = i + 1; j < MAXPOPS; j++) {
                /*printf("D%d%d = %1.2f, ",i,j,
                  0.5*KLDiv(i,j,P,NumAlleles)+0.5*KLDiv(j,i,P,NumAlleles) ); */
                div = 0.5 * KLDiv (i, j, P,LogP, NumAlleles, reps) + 0.5 * KLDiv (j, i, P,LogP,
                        NumAlleles, reps);
                fprintf (file, "%1.3f ", div);
                /*add extra spaces when i and j are large, to agree with spacing of banner */
                for (k = 2; k < ((int) log10 (i)) + ((int) log10 (j)); k++) {
                    fprintf (file, " ");
                }
            }
    } else
        /*print 2-D table */
    {
        fprintf (file,
                 "\nAllele-freq. divergence among pops (Kullback-Leibler distance)");
        if (reps > 1) {
            fprintf (file, ",\ncomputed using point estimates of P");
        }
        fprintf (file, ".\n\n");

        fprintf (file, "     ");
        for (j = 0; j < MAXPOPS; j++) {
            fprintf (file, "%2d    ", j + 1);
        }
        fprintf (file, "\n");

        for (i = 0; i < MAXPOPS; i++) {
            fprintf (file, "%2d  ", i + 1);
            for (j = 0; j < MAXPOPS; j++) {
                if (i == j) {
                    fprintf (file, "  -   ");
                } else {
                    fprintf (file, "%2.2f  ", KLDiv (i, j, P,LogP, NumAlleles, reps));
                }
            }
            fprintf (file, "\n");
        }
    }
}
/*---------------------------------------------------*/
void
UpdateSums (float *Q, float *QSum, int *Z, float *P, float *PSum,
            float *Fst, float *FstSum, int *NumAlleles,
            int *AncestDist, float *Epsilon, float *SumEpsilon,
            float *lambda, float *sumlambda, float *LocPrior,
            float *sumLocPrior, int LocPriorLen)
{
    int loc, ind, pop, allele, box, i;
      int line; 

    for (pop=0; pop<MAXPOPS; pop++) {
        sumlambda[pop] += lambda[pop];
    }

    for (ind = 0; ind < NUMINDS; ind++)
        for (pop = 0; pop < MAXPOPS; pop++) {
            QSum[QPos (ind, pop)] += Q[QPos (ind, pop)];
        }


    for (loc = 0; loc < NUMLOCI; loc++)
        for (pop = 0; pop < MAXPOPS; pop++)
            for (allele = 0; allele < NumAlleles[loc]; allele++) {
                PSum[PPos (loc, pop, allele)] += P[PPos (loc, pop, allele)];
            }

    if (FREQSCORR) {
        for (pop = 0; pop < MAXPOPS; pop++) {
            FstSum[pop] += Fst[pop];
        }

        for (loc = 0; loc < NUMLOCI; loc++)
            for (allele = 0; allele < NumAlleles[loc]; allele++) {
                SumEpsilon[EpsPos (loc, allele)] += Epsilon[EpsPos (loc, allele)];
            }

    }

    if (ANCESTDIST) {             /*store histogram of Q values for each individual */
        for (ind = 0; ind < NUMINDS; ind++)
            for (pop = 0; pop < MAXPOPS; pop++) {
                box = ((int) (Q[QPos (ind, pop)] * ((float) NUMBOXES)));
                /*printf("%1.3f__%d  ",Q[QPos(ind,pop)],box); */
                if (box == NUMBOXES) {
                    box = NUMBOXES - 1;    /*ie, Q = 1.000 */
                }
                AncestDist[AncestDistPos (ind, pop, box)]++;
            }
    }
    if (LOCPRIOR)
        for (i=0; i<LocPriorLen; i++) {
            sumLocPrior[i] += LocPrior[i];
        }

}


/*====================================================*/
/*---------------------------------------------------*/
void
PrintMainParams (FILE * file, int rep, int argc, char *argv[])
{
    int i;
    /*print values of most important parameters */

    fprintf (file, "\n");
    if (argc > 1) {
        fprintf (file, "Command line arguments:   ");
        for (i = 0; i < argc; i++) {
            fprintf (file, "%s ", argv[i]);
        }
        fprintf (file, "\n");
    }
    fprintf (file, "Input File:    %s\n", DATAFILE);

    if (file == stdout) {
        fprintf (file, "Output File:   %s_f\n", OUTFILE);
    }
    fprintf (file, "\n");
    fprintf (file, "Run parameters:\n");
    fprintf (file, "   %d individuals\n", NUMINDS);
    fprintf (file, "   %d loci\n", NUMLOCI);
    fprintf (file, "   %d populations assumed\n", MAXPOPS);
    fprintf (file, "   %d Burn-in period\n", BURNIN);
    fprintf (file, "   %d Reps\n", rep - BURNIN);

    if (USEPOPINFO) {
        fprintf (file, "USEPOPINFO turned on\n");
        fprintf (file, "MIGRPRIOR = %1.4f\n", MIGRPRIOR);
    }
    if (RECESSIVEALLELES) {
        fprintf(file,"RECESSIVE ALLELES model used\n");
    }
    if (NOADMIX) {
        fprintf (file, "NO ADMIXTURE model assumed\n");
    }
    if (STARTATPOPINFO) {
        fprintf (file, "STARTATPOPINFO turned on\n");
    }
    if (LOCPRIOR) {
        fprintf (file, "LOCPRIOR model used\n");
    }
    if (!(RANDOMIZE)) {
        fprintf (file, "RANDOMIZE turned off\n");
    }
    fprintf (file, "\n");

}
/*----------------------------------------------------*/
void
PrintAncestDist (FILE * file, int *AncestDist, int ind, int rep)
/*print credible region for each q in each population */
{
    int pop, box, low;
    float sum;
    /*print the whole histogram */
    /*fprintf(file,"\n");
      for (pop=0; pop<MAXPOPS; pop++)
      {
      fprintf(file,"%d: ",pop);
      for (box=0; box<NUMBOXES; box++)
      fprintf(file,"%d ",AncestDist[AncestDistPos(ind,pop,box)]);
      fprintf(file,"\n");
      } */

    fprintf (file, "    ");
    for (pop = 0; pop < MAXPOPS; pop++) {
        sum = 0;
        low = 0;
        for (box = 0; box < NUMBOXES; box++) {
            sum += AncestDist[AncestDistPos (ind, pop, box)];
            if ((low == 0) && (sum > (int) ((rep - BURNIN) * (1.0 - ANCESTPINT) / 2.0))) {
                /*printf("lo: sum = %d, value = %d\n",
                  (int) (rep-BURNIN)*(1.0-ANCESTPINT)/2.0); */
                fprintf (file, "(%1.3f,", (float) (box + 0.0) / NUMBOXES);
                low = 1;
            }
            if (sum > ((rep - BURNIN) * (1 + ANCESTPINT) / 2.0)) {
                /*{printf("hi: sum = %d, value = %d\n",
                  (int) (rep-BURNIN)*(1.0-ANCESTPINT)/2.0); */
                fprintf (file, "%1.3f) ", (float) (box + 1.0) / NUMBOXES);
                break;
            }
        }
    }
}
/*----------------------------------------------------*/
void
PrintGeogQ (FILE * file, int ind, int pop, float *UsePopProbs, int rep)
{
    int homepop, gen;

    homepop = pop - 1;
    fprintf (file, "%1.3f | ",
             (float) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN));
    for (pop = 0; pop < MAXPOPS; pop++)
        if (pop != homepop) {
            fprintf (file, "Pop %d: ", pop + 1);
            for (gen = 0; gen < GENSBACK + 1; gen++)
                fprintf (file, "%1.3f ",
                         (float) UsePopProbs[UsePPrPos (ind, pop, gen)] / (rep - BURNIN));
            fprintf (file, " | ");
        }

    fprintf (file, " ");
    if ((float) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN) < 0.1) {
        fprintf (file, "*");
    }
    if ((float) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN) < 0.3) {
        fprintf (file, "*");
    }
    if ((float) UsePopProbs[UsePPrPos (ind, homepop, 0)] / (rep - BURNIN) < 0.5) {
        fprintf (file, "*");
    }

}
/*----------------------------------------------------*/
void
PrintQ (FILE * file, int *Geno, int rep, float *QSum, struct IND *Individual,
        int *AncestDist, float *UsePopProbs,float *sumR)
/*print summary of Q to file */
{
    int ind, pop;
    int missing;                  /*(% missing data for this individual) */
    /*  int sum, low, box; */
    int MissingInd (int *Geno, int ind);


    fprintf (file, "\n");
    fprintf (file, "Inferred ancestry of individuals:\n");
    if (USEPOPINFO) {
        fprintf (file,
                 "Probability of being from assumed population | prob of other pops\n");
    }
    fprintf (file, "    ");
    if (LABEL) {
        fprintf (file, " %8s", "Label");
    }
    fprintf (file, " (%%Miss) ");
    if (POPDATA) {
        fprintf (file, "Pop");
    }
    if (!USEPOPINFO) {
        fprintf(file, ":  ");
    }
    if (LINKAGE && INDIVIDUALR) {
        fprintf(file, " Ind's r ");
    }
    if (!(USEPOPINFO)) {
        fprintf (file, "Inferred clusters");
        if (ANCESTDIST)
            fprintf (file, " (and %d%c probability intervals)",
                     (int) (100.0 * ANCESTPINT), '%');
    }
    fprintf (file, "\n");

    /*
      if (LABEL) fprintf(file,"    %8s ","");
      fprintf(file,"         ");
      if (POPDATA) fprintf(file,"   ");
      fprintf(file,":   ");
      for (pop=0; pop<MAXPOPS; pop++)
      fprintf(file,"%2d    ",pop+1);
      fprintf(file,"\n"); */

    for (ind = 0; ind < NUMINDS; ind++) {
        missing = (int) (100 * MissingInd (Geno, ind)) / (LINES * NUMLOCI);

        fprintf (file, "%3d ", ind + 1);
        if (LABEL) {
            fprintf (file, "%8s ", Individual[ind].Label);
        }
        if (missing < 10) {
            fprintf (file, " ");
        }
        fprintf (file, "  (%d) ", missing);
        fprintf (file, "  ");
        if (POPDATA) {
            fprintf (file, "%2d ", Individual[ind].Population);
        }
        fprintf (file, ":  ");
        if (LINKAGE && INDIVIDUALR) {
            fprintf(file, " %1.4f  ",(float)sumR[ind]/(rep-BURNIN));
        }
        if ((USEPOPINFO) && (Individual[ind].PopFlag)) {
            PrintGeogQ (file, ind, Individual[ind].Population, UsePopProbs, rep);
        }

        else {
            for (pop = 0; pop < MAXPOPS; pop++) {
                fprintf (file, "%1.3f ", (float) QSum[QPos (ind, pop)] / (rep - BURNIN));
            }

            if (ANCESTDIST) { /*Print the credible intervals for ancestry coeffs */
                PrintAncestDist (file, AncestDist, ind, rep);
            }
        }
        fprintf (file, "\n");
    }
}
/*-----------------------------------------------------*/
void
PrintQFile (int rep, float *QSum, struct IND *Individual, float *UsePopProbs)
/*Produce a file that contains only the label and popinfo (if any) plus Q-hat */
{
    char outname[STRLEN + 20];
    FILE *QHat;
    int ind, pop;
    /*float *QProbs; *//*[MAXPOPS] */

    sprintf (outname, "%s_q", OUTFILE);
    QHat = fopen (outname, "w");
    if (QHat == NULL) {
        printf ("WARNING: Unable to open output file %s.\n", outname);
    }

    /*QProbs = calloc(MAXPOPS,sizeof(float));
      if (QProbs==NULL) printf("Warning: unable to assign memory in PrintQFile\n"); */

    for (ind = 0; ind < NUMINDS; ind++) {
        if (LABEL) {
            fprintf (QHat, "%12s ", Individual[ind].Label);
        } else {
            fprintf (QHat, "%4d ", ind + 1);
        }
        if (POPDATA) {
            fprintf (QHat, "%2d ", Individual[ind].Population);
        }
        /*if ((USEPOPINFO)&&(Individual[ind].PopFlag))
          QFromUsePop(ind,QProbs,UsePopProbs,Individual[ind].Population-1);
          else
          for (pop=0; pop<MAXPOPS; pop++)
          QProbs = (float) QSum[QPos(ind,pop)]/(rep-BURNIN); */
        for (pop = 0; pop < MAXPOPS; pop++) {
            fprintf (QHat, "%1.4f ", (float) QSum[QPos (ind, pop)] / (rep - BURNIN));
        }
        fprintf (QHat, "\n");
    }
    fclose (QHat);
    /*free(QProbs); */
}
/*-----------------------------------------------------*/
void
PrintP (FILE * file, int rep, int *Geno, float *PSum, int *Translation,
        int *NumAlleles,
        float *SumEpsilon, char *Markername)
/*print summary of P to file */
{
    int loc, pop, allele;
    int MissingLoc (int *Geno, int loc);


    fprintf (file, "\n\nEstimated Allele Frequencies in each cluster\n");
    if (FREQSCORR) {
        fprintf (file, "First column gives estimated ancestral frequencies\n");
    }
    fprintf (file, "\n\n");

    for (loc = 0; loc < NUMLOCI; loc++) {
        fprintf (file, "Locus %d : ", loc + 1);
        if (MARKERNAMES) {
            PrintGeneName(file, loc, Markername);
        }
        fprintf (file, "\n");
        fprintf (file, "%d alleles\n", NumAlleles[loc]);
        fprintf (file, "%2.1f%c missing data\n",
                 (float) 100.0 * MissingLoc (Geno, loc) / (LINES * NUMINDS),
                 '%');  /* JKP changed 2 to LINES (4/06/09) */
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
            /*This fine piece of programming uses the allele 29697 to represent the
              recessive allele in the event that it is not observed in the data; used
              only when the RECESSIVEALLELES model is turned on. This is coded in datain.c*/
            if (RECESSIVEALLELES
                    && Translation[TransPos (loc, allele)] == 29697) {
                fprintf (file, "Null   ");
            }

            else {
                fprintf (file, "%4d   ", Translation[TransPos (loc, allele)]);
            }

            if (FREQSCORR)
                fprintf (file, "(%1.3f) ",
                         (float) SumEpsilon[EpsPos (loc, allele)] / (rep - BURNIN));

            for (pop = 0; pop < MAXPOPS; pop++) {
                fprintf (file, "%1.3f ", (float) PSum[PPos (loc, pop,
                                                       allele)] / (rep - BURNIN));
            }
            fprintf (file, "\n");
        }
        fprintf (file, "\n");
    }
}

int get_location_num(int loc, struct IND *individual)
{
    int ind;
    for (ind=0; ind<NUMINDS; ind++)
        if (individual[ind].myloc==loc) {
            return individual[ind].Location;
        }
    printf("error in get_location_num: can't find location %i\n", loc);
    exit(-1);
    return -1;
}

/*-----------------------------------------------------*/
void
PrintSums (FILE * file, int rep, float sumlikes,
           float sumsqlikes, float *FstSum, float *sumAlpha,
           float *sumlambda, float *sumR,float *varR,
           struct IND *Individual, float *sumLocPrior, int LocPriorLen,
           float DIC)
/*print current value of some averages to file */
{
    int ind, pop, locnum, loc;
    float sumrecs = 0.0;

    fprintf (file, "--------------------------------------------\n");
    if (COMPUTEPROB) {
        if (rep - BURNIN > 2) {

            fprintf (file, "Estimated Ln Prob of Data   = %1.1f\n",
                     EstLogProb (sumlikes, sumsqlikes, rep - BURNIN));
            fprintf (file, "Mean value of ln likelihood = %1.1f\n",
                     (float) sumlikes / (rep - BURNIN));
            fprintf (file, "Variance of ln likelihood   = %1.1f\n",
                     SampleVar (sumsqlikes, sumlikes, (rep - BURNIN)));
            /*          fprintf (file, "DIC = %.3f\n", DIC); */
        }

        else
            fprintf (file, "Mean value of ln likelihood = %1.1f\n",
                     (float) sumlikes / (rep - BURNIN));     /*print this line in either case */
    }

    if (((!(NOADMIX)) && (!(NOALPHA)) && (MAXPOPS > 1))) {
        if (POPALPHAS) {
            fprintf (file, "\n");
            for (pop = 0; pop < MAXPOPS; pop++)
                fprintf (file, "Mean value of alpha_%d       = %1.4f\n", pop + 1,
                         (float) sumAlpha[pop] / (rep - BURNIN));
        } else {
            fprintf (file, "Mean value of alpha         = %1.4f\n",
                     (float) sumAlpha[0] / (rep - BURNIN));
        }

        if (LOCPRIOR) {
            fprintf(file, "\nMean value of alpha_local for each location:\n");
            for (loc=0; loc<NUMLOCATIONS; loc++) {
                locnum = get_location_num(loc, Individual);
                fprintf(file, "\tlocation %2i:", locnum);
                for (pop=0; pop<MAXPOPS; pop++) {
                    fprintf(file, "  %1.4f", sumAlpha[AlphaPos(loc, pop)]/(rep - BURNIN));
                }
                fprintf(file, "\n");
            }
        }
    }

    if (INFERLAMBDA) {
        if (POPSPECIFICLAMBDA)
            for (pop=0; pop<MAXPOPS; pop++)
                fprintf (file, "\nMean value of lambda%d       = %1.4f\n",pop+1,
                         (float) sumlambda[pop] / (rep - BURNIN));
        else

            fprintf (file, "\nMean value of lambda        = %1.4f\n",
                     (float) sumlambda[0] / (rep - BURNIN));
    }
    if (LINKAGE) {
        if (!INDIVIDUALR) {
            fprintf (file, "Mean value of r              = %1.4f\n",
                     (float) sumR[0] / (rep - BURNIN));
            fprintf (file, "Standard deviation of r    = %1.4f\n",
                     sqrt((float) varR[0] / (float) (rep - BURNIN) - sumR[0] * sumR[0] /
                          (float) ((rep - BURNIN) * (rep - BURNIN))));
        } else {
            for (ind = 0; ind < NUMINDS; ind++) {
                sumrecs += sumR[ind];
            }
            fprintf (file, "Mean value of r           = %1.6f\n",
                     (float) sumrecs / NUMINDS / (rep - BURNIN));

        }
    }

    if (FREQSCORR) {
        if (ONEFST)
            fprintf (file, "Mean value of Fst           = %1.4f\n",
                     (float) FstSum[0] / (rep - BURNIN));
        else {
            fprintf (file, "\n");
            for (pop = 0; pop < MAXPOPS; pop++)
                fprintf (file, "Mean value of Fst_%d         = %1.4f\n", pop + 1,
                         (float) FstSum[pop] / (rep - BURNIN));
        }
    } else {
        fprintf (file, "Allele frequencies uncorrelated\n");
    }

    if (PFROMPOPFLAGONLY) {
        fprintf (file,
                 "Allele frequencies updated using individuals with POPFLAG=1 ONLY.\n");
    }

    fprintf (file, "\n");

    if (LOCPRIOR) {
        fprintf(file, "Mean value of r = %1.4f\n", sumLocPrior[0]/(rep-BURNIN));
        if (NOADMIX) {
            fprintf(file, "Mean value of nu = ");
            for (pop=0; pop<MAXPOPS; pop++) {
                fprintf(file, " %1.4f", sumLocPrior[LocPriorPos(NUMLOCATIONS,
                                                    pop)]/(rep-BURNIN));
            }
            fprintf(file, "\n");
            fprintf(file, "Mean value of gamma for each location:\n");
            for (loc=0; loc<NUMLOCATIONS; loc++) {
                locnum = get_location_num(loc, Individual);
                fprintf(file, "\tlocation %2i:", locnum);
                for (pop=0; pop<MAXPOPS; pop++) {
                    fprintf(file, "  %1.4f", sumLocPrior[LocPriorPos(loc, pop)]/(rep-BURNIN));
                }
                fprintf(file, "\n");
            }
        }
    }
}


/*-----------------------------------------------------*/
int
MissingLoc (int *Geno, int loc)
/*return the number of missing alleles at a locus */
{
    int ind, line;
    int sofar = 0;

    for (ind = 0; ind < NUMINDS; ind++)
        for (line = 0; line < LINES; line++) {
            if (Geno[GenPos (ind, line, loc)] == MISSING) {
                sofar++;
            }
        }

    return sofar;
}
/*-----------------------------------------------------*/
int
MissingInd (int *Geno, int ind)
/*return the number of missing alleles in an individual */
{
    int loc, line;
    int sofar = 0;

    for (loc = 0; loc < NUMLOCI; loc++)
        for (line = 0; line < LINES; line++) {
            if (Geno[GenPos (ind, line, loc)] == MISSING) {
                sofar++;
            }
        }

    return sofar;
}
/*----------------------------------------------------*/
void PrintGeneName(FILE * file, int loc, char *Markername)
{
    int i;
    for (i=0; i<GENELEN; i++) {
        if (Markername[MarkernamePos(loc,i)] != '\0') {
            fprintf(file,"%c",Markername[MarkernamePos(loc,i)]);
        } else {
            if (i==0) {
                fprintf(file,"XXX");
            }
            fprintf(file," ");
            break;
        }
    }
}
/*----------------------------------------------------*/
int EqualGeneNames(int loc1,int loc2,char *Markername)
/*returns 1 if the gene names are the same, otherwise 0*/
{
    int i;

    for (i=0; i<GENELEN; i++) {
        if (Markername[MarkernamePos(loc1,i)] != Markername[MarkernamePos(loc2,i)]) {
            return 0;
        }
        if (Markername[MarkernamePos(loc1,i)] == '\0') {
            return 1;
        }
    }
    return 1;
}
/*----------------------------------------------------*/
void
PrintMembership (FILE * file, float *QSum, struct IND *Individual)
/*Print summary of relationship between given populations,
  and cluster populations.  Requires POPDATA. An earlier, more
  complicated version of this is stored in backup.c */
{
    float *sumvals;              /*this array stores the sum of QSum for each of the
                                  the cluster populations 0..MAXPOPS-1, for individuals
                                  who are designated as coming from particular population. */
    float rowsum;                /*sum of the values in the array sumvals */
    int ind;
    int minpop;                   /*value of smallest (largest) population number */
    int maxpop;
    int pop, givenpop;
    int numfrompop;               /*number of individuals from each given pop */

    sumvals = calloc (MAXPOPS, sizeof (float));
    if (sumvals == NULL) {
        fprintf (file, "Error assigning memory in function PrintMembership\n");
    } else {

        if (POPDATA) {
            minpop = Individual[0].Population;        /*figure out min and max population names */
            maxpop = Individual[0].Population;
            for (ind = 1; ind < NUMINDS; ind++) {
                if (Individual[ind].Population < minpop) {
                    minpop = Individual[ind].Population;
                }
                if (Individual[ind].Population > maxpop) {
                    maxpop = Individual[ind].Population;
                }
            }

            if (sumvals == NULL) {
                fprintf (file, "Error assigning memory in function PrintMembership\n");
            } else {
                fprintf (file, "\n--------------------------------------------\n");
                fprintf (file, "Proportion of membership of each pre-defined\n");
                fprintf (file, " population in each of the %d clusters\n\n", MAXPOPS);

                fprintf (file, "Given    Inferred Clusters");
                for (pop = 3; pop < MAXPOPS; pop++) {
                    fprintf (file, "       ");
                }
                fprintf (file, "       Number of\n");

                fprintf (file, " Pop    ");
                for (pop = 0; pop < MAXPOPS; pop++) {
                    fprintf (file, "  %2d   ", pop + 1);
                }
                for (pop = MAXPOPS; pop < 3; pop++) {
                    fprintf (file, "       ");
                }
                fprintf (file, "   Individuals\n\n");

                for (givenpop = minpop; givenpop <= maxpop; givenpop++) {
                    for (pop = 0; pop < MAXPOPS; pop++) {
                        sumvals[pop] = 0.0;
                    }

                    numfrompop = 0;
                    for (ind = 0; ind < NUMINDS; ind++) {
                        if (givenpop == Individual[ind].Population) {
                            numfrompop++;
                            for (pop = 0; pop < MAXPOPS; pop++) {
                                sumvals[pop] += QSum[QPos (ind, pop)];
                            }
                        }
                    }
                    rowsum = 0.0;
                    for (pop = 0; pop < MAXPOPS; pop++) {
                        rowsum += sumvals[pop];
                    }

                    if (rowsum > 0.0) {
                        fprintf (file, "%3d:     ", givenpop);
                        for (pop = 0; pop < MAXPOPS; pop++) {
                            fprintf (file, "%1.3f  ", sumvals[pop] / rowsum);
                        }
                        for (pop = MAXPOPS; pop < 3; pop++) { /*number of individuals */
                            fprintf (file, "       ");
                        }
                        fprintf (file, "    %3d\n", numfrompop);
                    }

                }
            }
        } else
            /* no popdata */
        {
            for (pop = 0; pop < MAXPOPS; pop++) {
                sumvals[pop] = 0.0;
            }

            for (ind = 0; ind < NUMINDS; ind++)
                for (pop = 0; pop < MAXPOPS; pop++) {
                    sumvals[pop] += QSum[QPos (ind, pop)];
                }

            rowsum = 0.0;
            for (pop = 0; pop < MAXPOPS; pop++) {
                rowsum += sumvals[pop];
            }

            fprintf (file, "\n--------------------------------------------\n");
            fprintf (file, "Overall proportion of membership of the\n");
            fprintf (file, "sample in each of the %d clusters\n\n", MAXPOPS);

            fprintf (file, "Inferred Clusters\n");
            for (pop = 0; pop < MAXPOPS; pop++) {
                fprintf (file, " %2d    ", pop + 1);
            }
            fprintf (file, "\n");
            for (pop = 0; pop < MAXPOPS; pop++) {
                fprintf (file, "%1.3f  ", sumvals[pop] / rowsum);
            }
            fprintf (file, "\n\n");


        }

        free (sumvals);
        fprintf (file, "--------------------------------------------\n");
    }
}

float CalcDIC(int rep, float sumlikes, float* sumindlikes,
               float* indlikes_norm)
{
    float sumind=0.0, dic;
    int ind;
    for (ind=0; ind<NUMINDS; ind++) {
        sumind += log(sumindlikes[ind]/(rep-BURNIN))+indlikes_norm[ind];
    }
    dic = -4.0*sumlikes/(rep-BURNIN)+2.0*sumind;
    /*  return -4.0*sumlikes/(rep-BURNIN) + 2.0*log(sumindlikes)/(rep-BURNIN); */
    return dic;
}


/*====================================================*/
/*Melissa modified 7/12/07 to incorporate poppriors and DIC*/
void
OutPutResults (int *Geno, int rep, int savefreq,
               struct IND *Individual,
               float *PSum, float *QSum,
               float *FstSum, int *AncestDist, float *UsePopProbs,
               float sumlikes, float sumsqlikes, float *sumAlpha,
               float *sumR, float *varR,
               int *NumAlleles, int *Translation, int final,
               char *Markername, float *R, float *SumEpsilon,
               float *lambda, float *sumlambda,
               float *sumLocPrior, int LocPriorLen,
               float *sumindlikes, float *indlikes_norm,
               int argc, char *argv[])

/*final indicates that the program is terminating.  Stuff gets
  printed to the screen at this stage. */
{
    /*print a bunch of stuff to file, and a subset of this to the screen:
      P, Q, Fsts, Net distances, likelihood results, all parameters, amount of
      missing data for each individual... */

    char outname[STRLEN + 20];

    FILE *RESULTS;
    /*  int outputoption; */
    float DIC = 0.0;
    /*  DIC = CalcDIC(rep, sumlikes, sumindlikes, indlikes_norm); */
    if (final) {
        sprintf (outname, "%s_f", OUTFILE);
    } else {
        sprintf (outname, "%s_%d", OUTFILE, (rep - BURNIN) / savefreq);
    }

    RESULTS = fopen (outname, "w");
    if (RESULTS == NULL) {
        printf ("WARNING: Unable to open output file %s.\n", outname);
    } else {
        Welcome (RESULTS);
        if (final) {
            printf ("\nMCMC completed\n");
        }
        PrintMainParams (RESULTS, rep, argc, argv);
        PrintMembership (RESULTS, QSum, Individual);
        PrintNET (RESULTS, PSum, NumAlleles, rep - BURNIN, 1);
        /*if (final) PrintNET(stdout,PSum,NumAlleles,rep-BURNIN,1); */
        PrintSums (RESULTS, rep, sumlikes, sumsqlikes, FstSum, sumAlpha, sumlambda,
                   sumR,varR, Individual, sumLocPrior, LocPriorLen, DIC);

        PrintQ (RESULTS, Geno, rep, QSum, Individual, AncestDist, UsePopProbs,sumR);
        PrintP (RESULTS, rep, Geno, PSum, Translation, NumAlleles, SumEpsilon,
                Markername);
        if (final) {
            PrintQ (stdout, Geno, rep, QSum, Individual, AncestDist, UsePopProbs,sumR);
            if (PRINTQHAT) {
                PrintQFile (rep, QSum, Individual, UsePopProbs);
            }
            PrintMainParams (stdout, rep, argc, argv);
            PrintSums (stdout, rep, sumlikes, sumsqlikes, FstSum, sumAlpha, sumlambda,
                       sumR,varR, Individual, sumLocPrior, LocPriorLen, DIC);
            PrintMembership (stdout, QSum, Individual);
        }

        PrintAllParams (RESULTS);
        if (final) {
            printf ("Final results printed to file %s\n\n", outname);
        }
    }
    fclose (RESULTS);

}
