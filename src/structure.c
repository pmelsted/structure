/*=======================================================
        createCLBuffer(clDict,POPFLAGCL,sizeof(int)*NUMINDS,CL_MEM_READ_WRITE);

  STRUCTURE.C

  Program for inferring population structure using multilocus
  genotype data.

  Code written by Daniel Falush, Melissa Hubisz, and Jonathan Pritchard

  See additional details in README file.

  =========================================================*/

#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "ran.h"
#include "params.h"
#include "datain.h"
#include "output.h"
#include "Kernels.h"

/*Here come the updates */
#include "Updates/UpdateQ.h"
#include "Updates/UpdateZ.h"
#include "Updates/UpdateP.h"
#include "Updates/UpdateAlpha.h"
#include "Updates/UpdateGeno.h"
#include "Updates/UpdateLocPrior.h"
#include "Updates/UpdateEpsilon.h"
#include "Updates/UpdateFst.h"
#include "Updates/UpdateLambda.h"


#include "Init.h"
#include "Util.h"

/* Declared global, for the signal handling */
CLDict *clDict = NULL;

static void catch_function(int signo){
    handleCLErr(signo,clDict,"Signal caught, exiting!\n");
    exit(EXIT_FAILURE);
}

void initQ(double *Q){
    int ind, pop;
    for(ind=0;ind<NUMINDS;++ind){
        for (pop = 0; pop < MAXPOPS; pop++) {
            Q[QPos (ind, pop)] = (double) 1 / MAXPOPS;
        }
    }
}

int compareArrs(int Z[], int OldZ[], int total, char *names)
{
    int i;
    int notSame=0;
    int same =total;
    for(i = 0; i < total; i++) {
        if(Z[i] != OldZ[i]) {
            notSame++;
            same--;
        }
    }
    if(notSame > 0) {
        printf("%s differ!\n",names);
        printf("Difference : same: %d, not same %d\n",same,notSame);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
int compareDoubleArrs(double Z[], double OldZ[], int total, char *names)
{
    int i;
    int notSame=0;
    int same =total;
    for(i = 0; i < total; i++) {
        if(fabs(Z[i] - OldZ[i]) > 10e-8) {
            printf("%f,%f\n",Z[i],OldZ[i]);
            notSame++;
            same--;
        }
    }
    if(notSame > 0) {
        printf("%s differ!\n",names);
        printf("Difference : same: %d, not same %d\n",same,notSame);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

void compareZCLandZ(CLDict *clDict,int *OrigZ, double *Q, double *P,int *Geno,
                    double *randomArr)
{
    int *OldZ;
    int *Z;
    int ret;

    /*Allocate memory and copy the original over, so that it doesn't interfere */
    OldZ = calloc(ZSIZE,sizeof(int));
    Z = calloc(ZSIZE,sizeof(int));
    memcpy(Z,OrigZ,ZSIZE*sizeof(int));

    /*UpdateZ only writes to Z */
    UpdateZ(Z,Q,P,Geno,randomArr);
    memcpy(OldZ,Z,ZSIZE*sizeof(int));
    UpdateZCL(clDict,Z,Q,P,Geno,randomArr);
    ret = compareArrs(Z,OldZ,ZSIZE,"Z and old Z");
    if (ret == EXIT_FAILURE) {
        ReleaseCLDict(clDict);
        exit(EXIT_FAILURE);
    }

    free(OldZ);
    free(Z);
}


void comparePCLandP(CLDict *clDict,double *OrigP,
                    double *Epsilon, double *Fst, int *NumAlleles, int *Geno, int *Z,
                    double *lambda, struct IND *Individual, double *randomArr)
{

    double *OldP;
    double *P;
    int ret;

    OldP = calloc(PSIZE,sizeof(double));
    P = calloc(PSIZE,sizeof(double));
    memcpy(P,OrigP,PSIZE*sizeof(double));

    UpdateP (P, Epsilon, Fst, NumAlleles, Geno, Z, lambda, Individual,
             randomArr);
    memcpy(OldP,P,PSIZE*sizeof(double));
    UpdatePCL (clDict,P,Epsilon, Fst, NumAlleles, Geno, Z, lambda,
               Individual,
               randomArr);
    ret = compareDoubleArrs(P,OldP,PSIZE,"P and old P");
    if (ret == EXIT_FAILURE) {
        ReleaseCLDict(clDict);
        exit(EXIT_FAILURE);
    }
    free(OldP);
    free(P);
}

void initRandGens(CLDict *clDict){
    size_t global[1];
    int seed = rand();
    global[0] = NUMRANDGENS;
    printf("Seed: %d\n", seed);
    printf("Structure seed: %d\n", SEED);
    setKernelArgExplicit(clDict,InitRandGenKernel,sizeof(int),&seed,1);
    runKernel(clDict,InitRandGenKernel,1,global,"InitRandGen");
    finishCommands(clDict,"init rand gens");
}


/*=============MAIN======================================*/

int main (int argc, char *argv[])
{
    /*data--------- */
    int *Geno;                    /*NUMINDSxLINES: genotypes */
    double *R;                    /*NUMINDS */
    double *Mapdistance;          /*NUMLOCI */
    double *Phase;                /*NUMLOCI*NUMINDS */
    int *Phasemodel=NULL;         /*NUMINDS */
    char *Markername;             /*GENELEN*NUMLOCI */

    struct IND *Individual;       /*NUMINDS: records for each individual */
    int *Translation;             /*NUMLOCIxMAXALLELES: value of each coded allele */
    int *NumAlleles;              /*NUMLOCI: number of alleles at each locus */

    /* only used for recessive or inbreeding models: */
    int *PreGeno=
        NULL;           /*NUMINDSxLINESxNUMLOCI; diploid genotype if recessive alleles */
    int *Recessive=
        NULL;         /*NUMLOCI recessive allele at each locus, or -1 if there is none */


    /*Basic parameters */
    int *Z;                       /*NUMINDSx2xNUMLOCI: Z=pop of origin for each allele */
    int *Z1;
    double *Q;                    /*NUMINDSxMAXPOPS:  Q=ancestry of individuals */
    double *P;                    /*NUMLOCIxMAXPOPSxMAXALLELES: P=population allele freqs */
    double *Epsilon;              /*NUMLOCIxMAXALLELES: Dirichlet parameter for allele
                                  frequencies. This is either LAMBDA (if uncorrelated), or
                                  ancestral allele freqs if they are correlated */
    double *Fst;          /*MAXPOPS: Factor multiplied by epsilon under the Fst model */
    double *Alpha;                /*MAXPOPS: Dirichlet parameter for degree of admixture.
                                  Start this at ALPHA, and possibly change
                                  (if INFERALPHA==1) */
    double *lambda;                /*Dirichlet prior parameter for allele frequencies;
                                   start this at LAMBDA, and update if INFERLAMBDA*/
    double *sumlambda;
    /*Summaries */
    int    *NumLociPop;           /*NUMINDSxMAXPOPS: Number of alleles from each pop (by ind) */
    double *PSum;                 /*NUMLOCIxMAXPOPSxMAXALLELES: sum of AlFreqs */
    double *QSum;                 /*NUMINDSxMAXPOPS:  sum of Ancestries */
    double *FstSum;               /*MAXPOPS:  Sum of Fst */
    double *SumEpsilon=
        NULL;      /*NUMLOCIxMAXALLELES: sum of ancestral allele freqs*/
    double *sumAlpha;              /*MAXPOPS*/
    double *sumR;                 /*NUMINDS */
    double *varR;                 /*NUMINDS */
    double recomblikelihood=0.0;
    double *like;                  /*current likelihood value */
    double *sumlikes;              /*sum of likelihood values */
    double *sumsqlikes;            /*sum of squared likelihoods */

    int *popflags; /*The populationflags of individuals*/
    unsigned int *randGens;


    /*Melissa added 7/12/07 for calculating DIC*/
    double *sumIndLikes, *indLikesNorm;

    int    *AncestDist=
        NULL;      /*NUMINDS*MAXPOPS*NUMBOXES histogram of Q values */
    double *UsePopProbs=
        NULL;     /*NUMINDS*MAXPOPS*(GENSBACK+1) This is used when the
                                  population info is used.  It stores the probability that an
                                  individual has each of a specified set of ancestry amounts */
    /*loop variables-------------- */
    int rep;                      /*MCMC iterations so far */
    int savefreq;                 /*frequency of saving to file */
    int ind;

    /*Melissa's new variables added 7/12/07 to use priors based on sampling location*/
    double *LocPrior=NULL, *sumLocPrior=NULL, LocPriorLen=0;

    /* ======================= GPU Structure ======================== */
    /*Dict to that keeps track of CL info */
    /*CLDict *clDict = NULL;*/
    double * randomArr; /* array of random numbers */
    int POPFLAGINDS = 0;
    double invsqrtnuminds;
    enum BUFFER buffers[5];
    char         *names[5];
    size_t        sizes[5];
    void         *dests[5];

    double  *reduceresult;
    int *Numafrompopscl;
    int *Numlocipopscl;

    if (signal(SIGINT, catch_function) == SIG_ERR) {
        fputs("An error occurred while setting a signal handler.\n", stderr);
        return EXIT_FAILURE;
    }

    clDict = malloc(sizeof (*clDict));
    sumlikes = calloc(1,sizeof(double));
    sumsqlikes = calloc(1,sizeof(double));
    like = calloc(1,sizeof(double));
    /*=====Code for getting started=============================*/

    Welcome (stdout);             /*welcome */
    GetParams (0,argc,argv);      /*read in parameter values */

    CheckParamCombinations();     /*check that some parameter combinations are valid*/

    Mapdistance = calloc (NUMLOCI, sizeof (double));
    Phase = calloc (NUMLOCI * NUMINDS, sizeof (double));


    if (LINES ==2 && PHASED ==0) {
        Phasemodel=calloc(NUMINDS,sizeof(int));
        for (ind=0; ind<NUMINDS; ind++) {
            if (MARKOVPHASE) {
                Phasemodel[ind]=0;
            } else {
                Phasemodel[ind]=1;
            }
        }
    }

    lambda=calloc(MAXPOPS, sizeof (double));
    sumlambda=calloc(MAXPOPS, sizeof (double));

    Markername = calloc (GENELEN*NUMLOCI, sizeof (char));
    Geno = calloc (LINES * NUMLOCI * NUMINDS, sizeof (int));
    if (RECESSIVEALLELES) {
        PreGeno = calloc (LINES * NUMLOCI * NUMINDS, sizeof (int));
        Recessive = calloc (NUMLOCI, sizeof (int));
        if (PreGeno == NULL || Recessive == NULL) {
            printf ("Error (3) in assigning memory\n");
            Kill ();
        }
    }

    Individual = calloc (NUMINDS, sizeof (struct IND));
    if (Geno == NULL || Individual == NULL || Mapdistance == NULL
            || Markername == NULL) {
        printf ("Error in assigning memory (not enough space?)\n");
        Kill ();
    }
    Randomize(RANDOMIZE, &SEED);

    /*read in data file */
    if (RECESSIVEALLELES) {
        ReadInputFile(PreGeno, Mapdistance, Markername, Individual, Phase, Recessive);
    } else {
        ReadInputFile (Geno, Mapdistance, Markername, Individual, Phase, Recessive);
    }

    if (RECESSIVEALLELES) {
        MAXALLELES = FindMaxAlleles (PreGeno, Recessive);
    } else {
        MAXALLELES = FindMaxAlleles (Geno, Recessive);
    }


    /*=============set aside memory space=====================*/
    Translation = calloc (NUMLOCI * MAXALLELES, sizeof (int));
    NumAlleles = calloc (NUMLOCI, sizeof (int));
    Z = calloc (NUMINDS * LINES * NUMLOCI, sizeof (int));
    Z1 = calloc (NUMINDS * LINES * NUMLOCI, sizeof (int));
    Q = calloc (NUMINDS * MAXPOPS, sizeof (double));
    P = calloc (NUMLOCI * MAXPOPS * MAXALLELES, sizeof (double));
    R = calloc (NUMINDS, sizeof (double));
    sumR = calloc (NUMINDS, sizeof (double));
    varR = calloc (NUMINDS, sizeof (double));
    Epsilon = calloc (NUMLOCI * MAXALLELES, sizeof (double));
    if (FREQSCORR) {
        SumEpsilon = calloc (NUMLOCI * MAXALLELES, sizeof (double));
    }
    Fst = calloc (MAXPOPS, sizeof (double));
    FstSum = calloc (MAXPOPS, sizeof (double));
    NumLociPop = calloc (NUMINDS * MAXPOPS, sizeof (int));
    PSum = calloc (NUMLOCI * MAXPOPS * MAXALLELES, sizeof (double));
    QSum = calloc (NUMINDS * MAXPOPS, sizeof (double));


    if (ANCESTDIST) {
        AncestDist = calloc (NUMINDS * MAXPOPS * NUMBOXES, sizeof (int));
    }
    if (USEPOPINFO) {
        UsePopProbs = calloc (NUMINDS * MAXPOPS * (GENSBACK + 1), sizeof (double));
    }

    /*Melissa added 7/12/07*/
    if (LOCDATA>0 || LOCISPOP) {
        GetNumLocations(Individual);
    }

    /*Allocate the LocPrior vector.
      For no-admixture, it contains r, and the vectors nu and gamma.
      For admixture, it contains gamma.  The alphas_locals are stored with alpha global*/
    if (LOCPRIOR) {
        if (NOADMIX) {
            LocPriorLen = 1+MAXPOPS*(NUMLOCATIONS+1);
        } else {
            LocPriorLen=1;
        }
        LocPrior = malloc(LocPriorLen*sizeof(double));
        sumLocPrior = malloc(LocPriorLen*sizeof(double));
    }

    if (LOCPRIOR && NOADMIX==0) {
        Alpha = malloc(MAXPOPS*(NUMLOCATIONS+1)*sizeof(double));
        sumAlpha = malloc(MAXPOPS*(NUMLOCATIONS+1)*sizeof(double));
    } else {
        Alpha = calloc(MAXPOPS, sizeof (double));
        sumAlpha = calloc(MAXPOPS, sizeof (double));
    }

    /* this is for DIC */
    sumIndLikes = malloc(NUMINDS*sizeof(double));
    indLikesNorm = malloc(NUMINDS*sizeof(double));

    if ((Translation == NULL) || (NumAlleles == NULL) || (Z == NULL)
            || (Z1 == NULL) || (Q == NULL) ||
            (P == NULL) || (R == NULL) || (sumR == NULL) || (varR == NULL)
            || (Epsilon == NULL) ||
            (Fst == NULL) || (NumLociPop == NULL) ||
            (PSum == NULL) || (QSum == NULL) ||  (FstSum == NULL) ||
            ((ANCESTDIST) && (AncestDist == NULL)) ||
            ((USEPOPINFO) && (UsePopProbs == NULL))||(Alpha == NULL)||(sumAlpha==NULL)||
            ((FREQSCORR) && (SumEpsilon == NULL)) ||
            (LocPriorLen>0 && (LocPrior==NULL || sumLocPrior==NULL)) ||
            sumIndLikes==NULL || indLikesNorm==NULL) {

        printf ("Error in assigning memory (not enough space?)\n");
        FreeAll(Mapdistance, Phase, Phasemodel, lambda, sumlambda, Markername, Geno,
                PreGeno, Recessive,
                Individual, Translation, NumAlleles, Z, Z1, Q, P,  R, sumR, varR, Epsilon,
                SumEpsilon,
                Fst, FstSum, NumLociPop, PSum, QSum,  AncestDist, UsePopProbs, LocPrior,
                sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm, clDict);
        Kill ();
    }
    /*=========done setting aside memory space=====================*/

    /*initialize variables and arrays */
    Initialization (Geno, PreGeno, Individual, Translation, NumAlleles, Z, Z1,
                    Epsilon, SumEpsilon,
                    Fst, PSum, Q, QSum, FstSum, AncestDist, UsePopProbs, Alpha,
                    sumAlpha, sumR, varR, sumlikes, sumsqlikes, &savefreq, R, lambda,
                    sumlambda,Phase,Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes,
                    indLikesNorm, clDict);

    /* ==================== GPU Structure ==================== */

    /*Allocate an array of random numbers. Primarily used so that we can compare
     CL implementation to the original */
    randomArr = calloc(RANDSIZE,sizeof(double));
    randGens = calloc(NUMRANDGENS,sizeof(unsigned int));
    Numafrompopscl = calloc(NUMLOCI*MAXPOPS*MAXALLELES,sizeof(int));
    Numlocipopscl = calloc(NUMINDS*MAXPOPS,sizeof(int));

    /* ====== OpenCL initialized ====== */

    printf ("\n\n--------------------------------------\n\n");
    printf ("Finished initialization; starting MCMC \n");
    printf ("%d iterations + %d burnin\n\n", NUMREPS, BURNIN);

    /*=====Main MCMC loop=======================================*/

    writeBuffer(clDict,Numafrompopscl,sizeof(int) * NUMLOCI*MAXPOPS*MAXALLELES,NUMAFROMPOPSCL,"NumAFromPops");
    writeBuffer(clDict,Numlocipopscl,sizeof(int)*NUMINDS*MAXPOPS,NUMLOCIPOPSCL,"NUMLOCIPOPS");
    /* init buffers on GPU */
    writeBuffer(clDict,P,sizeof(double) * PSIZE,PCL,"P");
    writeBuffer(clDict,Z,sizeof(int)*ZSIZE,ZCL,"Z");
    writeBuffer(clDict,NumAlleles,sizeof(int) * NUMLOCI,NUMALLELESCL,"NumAlleles");
    writeBuffer(clDict,lambda,sizeof(double) * MAXPOPS,LAMBDACL,"LAMBDA");
    if(!RECESSIVEALLELES){
        writeBuffer(clDict,Geno,sizeof(int)*GENOSIZE,GENOCL,"Geno");
    }

    popflags = calloc(NUMINDS,sizeof(int));
    for(ind = 0; ind < NUMINDS;ind++){
        popflags[ind] = Individual[ind].PopFlag;
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            POPFLAGINDS++;
        }
    }

    printf("Setting updatealpha arg\n");
    setKernelArgExplicit(clDict,UpdateAlphaKernel,sizeof(int),&POPFLAGINDS,7);

    printf("Setting invsqrtnuminds arg\n");
    invsqrtnuminds = 1.0/sqrt(NUMINDS);
    setKernelArgExplicit(clDict,NonIndUpdateEpsilonKernel,sizeof(double),&invsqrtnuminds,6);

    writeBuffer(clDict,popflags,sizeof(int)*NUMINDS,POPFLAGCL,"popflags");

    writeBuffer(clDict,Alpha,sizeof(double) *MAXPOPS,ALPHACL,"Alpha");
    reduceresult = calloc(NUMINDS*NUMLOCI*MAXGROUPS,sizeof(double));
    if(reduceresult == NULL){
        printf("Failed to allocate reduce result\n");
    }
    printf("Writing reduce results\n");
    writeBuffer(clDict,reduceresult,sizeof(double)*MAXGROUPS*NUMINDS*NUMLOCI,REDUCERESULTSCL,"result");
    printf("Done writing reduce results\n");

    writeBuffer(clDict,sumAlpha, sizeof(double) * MAXPOPS,ALPHASUMCL, "alphasum");
    writeBuffer(clDict,sumlambda, sizeof(double) * MAXPOPS,LAMBDASUMCL, "lambdasum");

    writeBuffer(clDict,like, sizeof(double),LIKECL, "like");
    writeBuffer(clDict,sumsqlikes, sizeof(double),SUMSQLIKESCL, "sumsqlikes");
    writeBuffer(clDict,sumlikes, sizeof(double),SUMLIKESCL, "sumlikes");

    writeBuffer(clDict,QSum, sizeof(double) * QSIZE,QSUMCL, "qsum");
    writeBuffer(clDict,PSum, sizeof(double) * PSIZE,PSUMCL, "psum");
    writeBuffer(clDict,SumEpsilon, sizeof(double) * NUMLOCI*MAXALLELES,EPSILONSUMCL, "epssum");
    if(ANCESTDIST){
        writeBuffer(clDict,AncestDist, sizeof(int)*NUMINDS*MAXPOPS*NUMBOXES,ANCESTDISTCL, "ancest dist");
    }

    if (FREQSCORR) {
        writeBuffer(clDict,Fst,sizeof(double) * MAXPOPS,FSTCL,"FST");
        writeBuffer(clDict,FstSum, sizeof(double) * MAXPOPS,FSTSUMCL, "fstsum");
        writeBuffer(clDict,Epsilon,sizeof(double) * NUMLOCI*MAXALLELES,EPSILONCL,
                    "EPSILON");
    }
    /*printf("%d, %d\n",INFERALPHA,INFERLAMBDA);*/
    /*printf("%d\n",USEPOPINFO);*/
    /*printf("%d\n",RECESSIVEALLELES);*/
    /*printf("%d\n",LINKAGE);*/
    /*printf("%d\n",COMPUTEPROB);*/
    /*printf("%d\n",POPALPHAS);*/
    /*printf("%d\n",NOADMIX);*/
    /*printf("%d\n",LOCPRIOR);*/
    /*handleCLErr(1,clDict,"heyhey");*/

    /*Initialize Q */
    initQ(Q);
    initRandGens(clDict);
    for (rep = 0; rep < (NUMREPS + BURNIN); rep++) {

        /*FillArrayWithRandomCL(clDict,randomArr,NUMLOCI*MAXALLELES*MAXPOPS*MAXRANDOM);*/
        /*FillArrayWithRandom(randomArr,NUMLOCI*MAXALLELES*MAXPOPS*MAXRANDOM);*/
        /*FillArrayWithRandom(randomArr,RANDSIZE);*/

        if(DEBUGCOMPARE) {
            readBuffer(clDict,randomArr,
                    sizeof(double) * NUMLOCI*MAXALLELES*MAXPOPS*MAXRANDOM,RANDCL,
                    "randomArr");
            comparePCLandP(clDict,P,Epsilon, Fst, NumAlleles, Geno, Z,
                           lambda, Individual, randomArr);
        }
        if (USEWORKINGCL) {
            /* clear buffer */
            writeBuffer(clDict,Numafrompopscl,sizeof(int) * NUMLOCI*MAXPOPS*MAXALLELES,NUMAFROMPOPSCL,"NumAFromPops");
            UpdatePCL (clDict,P, Epsilon, Fst, NumAlleles, Geno, Z, lambda,
                       Individual,
                       randomArr);
        }  else {
            readBuffer(clDict,randomArr,
                    sizeof(double) * NUMLOCI*MAXALLELES*MAXPOPS*MAXRANDOM,RANDCL,
                    "randomArr");
            UpdateP (P, Epsilon, Fst, NumAlleles, Geno, Z, lambda, Individual,
                     randomArr);
        }

        /* Update Q */
        /*FillArrayWithRandomCL(clDict,randomArr,NUMINDS+NUMINDS*MAXRANDOM);*/
        if (LINKAGE && rep >= ADMBURNIN) {
            UpdateQMetroRecombine (Geno, Q, Z, P, Alpha, rep,
                                   Individual, Mapdistance, R, Phase,Phasemodel,randomArr);
        } else {
            writeBuffer(clDict,Numlocipopscl,sizeof(int)*NUMINDS*MAXPOPS,NUMLOCIPOPSCL,"NUMLOCIPOPS");
            UpdateQCL (clDict,Geno, PreGeno, Q, P, Z, Alpha, rep, Individual, UsePopProbs,
                     Recessive, LocPrior,randomArr);
        }

        if (LOCPRIOR && UPDATELOCPRIOR) {
            UpdateLocPrior(Q, LocPrior, Alpha, Individual);
        }

        if (RECESSIVEALLELES) {
            UpdateGeno (PreGeno, Geno, P, Z, Recessive, NumAlleles, Q);
            writeBuffer(clDict,Geno,sizeof(int) * GENOSIZE,GENOCL,"Geno");
            /*The Zs are not correct after UpdateGeno, until UpdateZ is run */
        }

        if (LINKAGE && rep > ADMBURNIN) {
            if (!INDIVIDUALR) {
                recomblikelihood = UpdateZandSingleR(Z,  Q, P, Geno,
                                                     R, Mapdistance, rep, Phase, Z1,Phasemodel, rep+1 > BURNIN? sumIndLikes : NULL,
                                                     indLikesNorm);
            } else {
                recomblikelihood = UpdateZandR(Z,  Q, P, Geno, R,
                                               Mapdistance, rep, Phase, Z1,Phasemodel, rep+1 > BURNIN ? sumIndLikes:NULL,
                                               indLikesNorm);
            }
        } else {
            /*FillArrayWithRandomCL(clDict,randomArr,NUMINDS*NUMLOCI*LINES);*/
            if (DEBUGCOMPARE) {
                readBuffer(clDict,randomArr,
                        sizeof(double) * NUMINDS*NUMLOCI*LINES,RANDCL,
                        "randomArr");
                compareZCLandZ(clDict,Z,Q,P,Geno,randomArr);
            }
            if (USEWORKINGCL) {
                UpdateZCL (clDict,Z,  Q, P, Geno,randomArr);
                /* Not needed */
                /*readBuffer(clDict,Z,sizeof(int)*ZSIZE,ZCL,"Z");*/
            } else {
                readBuffer(clDict,randomArr,
                        sizeof(double) * NUMINDS*NUMLOCI*LINES,RANDCL,
                        "randomArr");
                UpdateZ (Z,  Q, P, Geno,randomArr);
            }
            /*      printf("done updatez alpha[2]=%e\n", Alpha[2]); */
        }


        if (LOCPRIOR && NOADMIX==0) {
            UpdateAlphaLocPrior(Q, Alpha, LocPrior, Individual);
        } else if (INFERALPHA) {

            UpdateAlphaCL (clDict,Q, Alpha, Individual, rep,POPFLAGINDS);

            /* readBuffer(clDict,Q, sizeof(double) * QSIZE,QCL, "Q"); */
            /* readBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha"); */

            /* UpdateAlpha(Q, Alpha, Individual, rep); */

            /* writeBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha"); */
        }

        if (INFERLAMBDA) {
            readBuffer(clDict,P,sizeof(double) * PSIZE,PCL,"P");
            if  (POPSPECIFICLAMBDA) {
                UpdatePopLambda(P,lambda,NumAlleles);
            } else {
                UpdateLambda (P,Epsilon,lambda, NumAlleles);
            }
            writeBuffer(clDict,lambda,sizeof(double) * MAXPOPS,LAMBDACL,"LAMBDA");
        }

        if (FREQSCORR) {
            UpdateEpsilonCL(clDict,P,Epsilon,Fst,NumAlleles,lambda[0]);

            UpdateFstCL (clDict,Epsilon, Fst, P, NumAlleles);

            /* readBuffer(clDict,P, sizeof(double) * PSIZE,PCL, "P"); */
            /* readBuffer(clDict,Fst,sizeof(double) * MAXPOPS,FSTCL,"FST"); */
            /* readBuffer(clDict,Epsilon,sizeof(double) * NUMLOCI*MAXALLELES,EPSILONCL,"eps"); */
            /* readBuffer(clDict,NumAlleles,sizeof(int) * NUMLOCI,NUMALLELESCL,"NumAlleles"); */

            /* UpdateEpsilon(P,Epsilon,Fst,NumAlleles,lambda[0]); */

            /* UpdateFst (Epsilon, Fst, P, NumAlleles); */

            /* writeBuffer(clDict,Epsilon,sizeof(double) * NUMLOCI*MAXALLELES,EPSILONCL,"eps"); */
            /* writeBuffer(clDict,Fst,sizeof(double) * MAXPOPS,FSTCL,"FST"); */
        }


        /*====book-keeping stuff======================*/
        if (rep + 1 > BURNIN) {
            /*buffers[0] = PCL;
            names[0] = "P"; dests[0] = P; sizes[0] = sizeof(double) * PSIZE;
            buffers[1] = QCL;
            names[1] = "Q"; dests[1] = Q; sizes[1] = sizeof(double) * QSIZE;

            readBuffers(clDict,dests,sizes,buffers,names,2);*/
            /*readBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha");
            if(rep % 100 == 0){
            printf("%f",Alpha[0]);
            }*/
            DataCollectionCL (clDict,Geno, PreGeno, Q, QSum, Z, Z1,  P, PSum,
                            Fst, FstSum, NumAlleles,
                            AncestDist, Alpha, sumAlpha, sumR, varR, like,
                            sumlikes, sumsqlikes, R, Epsilon,SumEpsilon,recomblikelihood,
                            lambda, sumlambda, Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes,
                            indLikesNorm, rep);
        }

        if ((savefreq) && ((rep + 1) > BURNIN)
                && (((rep + 1 - BURNIN) % savefreq) == 0)
                && ((rep + 1) != NUMREPS + BURNIN)) {
            readBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha");
            readBuffer(clDict,Fst, sizeof(double) * MAXPOPS,FSTCL, "fst");
            readBuffer(clDict,sumAlpha, sizeof(double) * MAXPOPS,ALPHASUMCL, "alphasum");
            readBuffer(clDict,sumlambda, sizeof(double) * MAXPOPS,LAMBDASUMCL, "lambdasum");
            readBuffer(clDict,FstSum, sizeof(double) * MAXPOPS,FSTSUMCL, "fstsum");
            readBuffer(clDict,QSum, sizeof(double) * QSIZE,QSUMCL, "Qsum");
            readBuffer(clDict,PSum, sizeof(double) * PSIZE,PSUMCL, "Psum");
            readBuffer(clDict,SumEpsilon, sizeof(double) * NUMLOCI*MAXALLELES,EPSILONSUMCL, "epssum");
            readBuffer(clDict,like, sizeof(double),LIKECL, "like");
            readBuffer(clDict,sumlikes, sizeof(double),SUMLIKESCL, "sumlike");
            readBuffer(clDict,sumsqlikes, sizeof(double),SUMSQLIKESCL, "sumsqlikes");
            OutPutResults (Geno, rep + 1, savefreq, Individual, PSum, QSum,
                           FstSum, AncestDist, UsePopProbs, *sumlikes,
                           *sumsqlikes, sumAlpha, sumR, varR,
                           NumAlleles, Translation, 0, Markername, R,
                           SumEpsilon,
                           lambda,sumlambda,sumLocPrior, LocPriorLen,
                           sumIndLikes, indLikesNorm, argc,argv);
        }


        if (PRINTLIKES) {
            readBuffer(clDict,like, sizeof(double),LIKECL, "like");
            PrintLike (*like, rep, Geno, PreGeno, Q, P,recomblikelihood);
        }

        if (((rep + 1) % UPDATEFREQ) == 0) {
            readBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha");
            readBuffer(clDict,Fst, sizeof(double) * MAXPOPS,FSTCL, "fst");
            readBuffer(clDict,like, sizeof(double),LIKECL, "like");
            readBuffer(clDict,sumlikes, sizeof(double),SUMLIKESCL, "sumlike");
            readBuffer(clDict,sumsqlikes, sizeof(double),SUMSQLIKESCL, "sumsqlikes");
            /*readBuffer(clDict,sumAlpha, sizeof(double) * MAXPOPS,ALPHASUMCL, "alphasum");
            readBuffer(clDict,sumlambda, sizeof(double) * MAXPOPS,LAMBDASUMCL, "lambdasum");
            readBuffer(clDict,FstSum, sizeof(double) * MAXPOPS,FSTSUMCL, "fstsum");
            readBuffer(clDict,QSum, sizeof(double) * QSIZE,QSUMCL, "Qsum");
            readBuffer(clDict,PSum, sizeof(double) * PSIZE,PSUMCL, "Psum");
            readBuffer(clDict,SumEpsilon, sizeof(double) * NUMLOCI*MAXALLELES,EPSILONSUMCL, "epssum");*/
            PrintUpdate (rep + 1, Geno, PreGeno, Alpha, Fst, P, Q, *like,
                         *sumlikes, *sumsqlikes, NumAlleles, R, lambda,Individual,
                         recomblikelihood, Recessive, LocPrior, LocPriorLen);
        }
    }

    /*====final book-keeping====================================*/
    if ((rep % UPDATEFREQ) != 0) {
        readBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha");
        readBuffer(clDict,Fst, sizeof(double) * MAXPOPS,FSTCL, "fst");
        readBuffer(clDict,like, sizeof(double),LIKECL, "like");
        readBuffer(clDict,sumsqlikes, sizeof(double),SUMSQLIKESCL, "sumsqlikes");
            readBuffer(clDict,sumlikes, sizeof(double),SUMLIKESCL, "sumlikes");
        /*readBuffer(clDict,sumAlpha, sizeof(double) * MAXPOPS,ALPHASUMCL, "alphasum");
        readBuffer(clDict,sumlambda, sizeof(double) * MAXPOPS,LAMBDASUMCL, "lambdasum");
        readBuffer(clDict,FstSum, sizeof(double) * MAXPOPS,FSTSUMCL, "fstsum");
        readBuffer(clDict,QSum, sizeof(double) * QSIZE,QSUMCL, "Qsum");
        readBuffer(clDict,PSum, sizeof(double) * PSIZE,PSUMCL, "Psum");
        readBuffer(clDict,SumEpsilon, sizeof(double) * NUMLOCI*MAXALLELES,EPSILONSUMCL, "epssum");*/
        PrintUpdate (rep, Geno, PreGeno, Alpha, Fst, P, Q, *like, *sumlikes,
                     *sumsqlikes, NumAlleles,R, lambda, Individual,recomblikelihood,
                     Recessive, LocPrior, LocPriorLen);
    }


    readBuffer(clDict,Alpha, sizeof(double) * MAXPOPS,ALPHACL, "alpha");
    readBuffer(clDict,Fst, sizeof(double) * MAXPOPS,FSTCL, "fst");
    readBuffer(clDict,sumAlpha, sizeof(double) * MAXPOPS,ALPHASUMCL, "alphasum");
    readBuffer(clDict,sumlambda, sizeof(double) * MAXPOPS,LAMBDASUMCL, "lambdasum");
    readBuffer(clDict,FstSum, sizeof(double) * MAXPOPS,FSTSUMCL, "fstsum");
    readBuffer(clDict,QSum, sizeof(double) * QSIZE,QSUMCL, "Qsum");
    readBuffer(clDict,PSum, sizeof(double) * PSIZE,PSUMCL, "Psum");
    readBuffer(clDict,SumEpsilon, sizeof(double) * NUMLOCI*MAXALLELES,EPSILONSUMCL, "epssum");
    readBuffer(clDict,like, sizeof(double),LIKECL, "like");
    readBuffer(clDict,sumsqlikes, sizeof(double),SUMSQLIKESCL, "sumsqlikes");
    readBuffer(clDict,sumlikes, sizeof(double),SUMLIKESCL, "sumlikes");
    OutPutResults (Geno, rep, savefreq, Individual, PSum, QSum,
                   FstSum, AncestDist, UsePopProbs,
                   *sumlikes, *sumsqlikes,
                   sumAlpha, sumR, varR, NumAlleles, Translation, 1,
                   Markername, R, SumEpsilon,
                   lambda,sumlambda,sumLocPrior, LocPriorLen,
                   sumIndLikes, indLikesNorm,
                   argc,argv);

    /*=====Closing everything down==============================*/
    FreeAll(Mapdistance, Phase, Phasemodel, lambda, sumlambda, Markername, Geno,
            PreGeno, Recessive,
            Individual, Translation, NumAlleles, Z, Z1, Q, P, R, sumR, varR, Epsilon,
            SumEpsilon,
            Fst, FstSum, NumLociPop, PSum, QSum,  AncestDist, UsePopProbs, LocPrior,
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm, clDict);
    free(randomArr);
    free(randGens);
    free(like);
    free(sumsqlikes);
    free(sumlikes);
    free(reduceresult);
    printf("Structure seed: %d\n", SEED);
    return (0);
}

/*==========================================

  Notes:


  ============================================*/
