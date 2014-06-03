/*=======================================================

  STRUCTURE.C

  Program for inferring population structure using multilocus
  genotype data.

  Code written by Daniel Falush, Melissa Hubisz, and Jonathan Pritchard

  See additional details in README file.

  =========================================================*/

#include <stdio.h>
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

void compareZAndOldZ(int Z[], int OldZ[]){
    int i;
    int notSame=0;
    int same =ZSIZE;
    for(i = 0; i < ZSIZE;i++){
        if(Z[i] != OldZ[i]){
            notSame++;
            same--;
        }
    }
    if(notSame > 0){
       printf("Z and Old Z differ!\n");
       printf("Difference between Z and OldZ: same: %d, not same %d\n",same,notSame);
       exit(EXIT_FAILURE);     
    }
}

void compareZCLandZ(CLDict *clDict,int *OrigZ, double *Q, double *P,int *Geno, double *randomArr){
    int *OldZ;
    int *Z;
    int i;
    /*Allocate memory and copy the original over, so that it doesn't interfere */
    OldZ = calloc(ZSIZE,sizeof(int));
    Z = calloc(ZSIZE,sizeof(int));
    memcpy(Z,OrigZ,ZSIZE*sizeof(int));
    
    /*UpdateZ only writes to Z */
    UpdateZ(Z,Q,P,Geno,randomArr);
    memcpy(OldZ,Z,ZSIZE*sizeof(int));
    UpdateZCL(clDict,Z,Q,P,Geno,randomArr);
    compareZAndOldZ(Z,OldZ);

    free(OldZ);
    free(Z);
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
  int *PreGeno=NULL;           /*NUMINDSxLINESxNUMLOCI; diploid genotype if recessive alleles */
  int *Recessive=NULL;         /*NUMLOCI recessive allele at each locus, or -1 if there is none */


  /*Basic parameters */
  int *Z;                       /*NUMINDSx2xNUMLOCI: Z=pop of origin for each allele */
  int *Z1;
  double *Q;                    /*NUMINDSxMAXPOPS:  Q=ancestry of individuals */
  double *P;                    /*NUMLOCIxMAXPOPSxMAXALLELES: P=population allele freqs */
  double *LogP;                 /*NUMLOCIxMAXPOPSxMAXALLELES: log of P, used to prevent underflow */
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
  double *SumEpsilon=NULL;      /*NUMLOCIxMAXALLELES: sum of ancestral allele freqs*/
  double *sumAlpha;              /*MAXPOPS*/
  double *sumR;                 /*NUMINDS */
  double *varR;                 /*NUMINDS */
  double recomblikelihood=0.0;
  double like;                  /*current likelihood value */
  double sumlikes;              /*sum of likelihood values */
  double sumsqlikes;            /*sum of squared likelihoods */


  /*Melissa added 7/12/07 for calculating DIC*/
  double *sumIndLikes, *indLikesNorm;

  int    *AncestDist=NULL;      /*NUMINDS*MAXPOPS*NUMBOXES histogram of Q values */
  double *UsePopProbs=NULL;     /*NUMINDS*MAXPOPS*(GENSBACK+1) This is used when the
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
  int ret;
  int i = 0;
  CLDict *clDict = NULL;
  char *names[6];
  char *vals[6];
  double * randomArr; /* array of random numbers */

  clDict = malloc(sizeof (*clDict));
  /*=====Code for getting started=============================*/

  Welcome (stdout);             /*welcome */
  GetParams (0,argc,argv);      /*read in parameter values */

  CheckParamCombinations();     /*check that some parameter combinations are valid*/
   

  Mapdistance = calloc (NUMLOCI, sizeof (double));
  Phase = calloc (NUMLOCI * NUMINDS, sizeof (double));


  if (LINES ==2 && PHASED ==0) {
    Phasemodel=calloc(NUMINDS,sizeof(int));
    for (ind=0;ind<NUMINDS;ind++) {
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
      printf ("Error (3) in assigning memory\n");Kill ();
    }
  }

  Individual = calloc (NUMINDS, sizeof (struct IND));
  if (Geno == NULL || Individual == NULL || Mapdistance == NULL || Markername == NULL) {
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
  LogP = calloc(NUMLOCI * MAXPOPS * MAXALLELES, sizeof(double));
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

  if ((Translation == NULL) || (NumAlleles == NULL) || (Z == NULL) || (Z1 == NULL) || (Q == NULL) ||
      (P == NULL) || (LogP==NULL) || (R == NULL) || (sumR == NULL) || (varR == NULL) || (Epsilon == NULL) ||
      (Fst == NULL) || (NumLociPop == NULL) ||
      (PSum == NULL) || (QSum == NULL) ||  (FstSum == NULL) ||
      ((ANCESTDIST) && (AncestDist == NULL)) ||
      ((USEPOPINFO) && (UsePopProbs == NULL))||(Alpha == NULL)||(sumAlpha==NULL)||
      ((FREQSCORR) && (SumEpsilon == NULL)) ||
      (LocPriorLen>0 && (LocPrior==NULL || sumLocPrior==NULL)) ||
      sumIndLikes==NULL || indLikesNorm==NULL) {

    printf ("Error in assigning memory (not enough space?)\n");
    FreeAll(Mapdistance, Phase, Phasemodel, lambda, sumlambda, Markername, Geno, PreGeno, Recessive,
            Individual, Translation, NumAlleles, Z, Z1, Q, P, LogP, R, sumR, varR, Epsilon, SumEpsilon,
            Fst, FstSum, NumLociPop, PSum, QSum,  AncestDist, UsePopProbs, LocPrior,
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm, clDict);
    Kill ();
  }
  /*=========done setting aside memory space=====================*/

  /*initialize variables and arrays */
  Initialization (Geno, PreGeno, Individual, Translation, NumAlleles, Z, Z1, Epsilon, SumEpsilon,
                  Fst, PSum, Q, QSum, FstSum, AncestDist, UsePopProbs, Alpha,
                  sumAlpha, sumR, varR, &sumlikes, &sumsqlikes, &savefreq, R, lambda,
                  sumlambda,Phase,Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes, indLikesNorm, clDict);

  
  /* ==================== GPU Structure ==================== */  
  /* compile OpenCL kernels */ 
  /*Define the constants in the kernels */

  for(i = 0; i < 6; ++i){
    vals[i] = calloc(255, sizeof(char));
  }

  names[0] = "%maxpops%"; sprintf(vals[0],"%d",MAXPOPS);
  names[1] = "%missing%"; sprintf(vals[1],"%d",MISSING);
  names[2] = "%maxalleles%"; sprintf(vals[2],"%d",MAXALLELES);
  names[3] = "%numloci%"; sprintf(vals[3],"%d",NUMLOCI);
  names[4] = "%lines%"; sprintf(vals[4],"%d",LINES);
  names[5] = "%numinds%"; sprintf(vals[5],"%d",NUMINDS);

  ret = CompileKernels(clDict,names,vals,6);
  for(i = 0; i < 6; ++i){
    free(vals[i]);
  }

  if(ret != EXIT_SUCCESS){
      printf("Kernels failed to compile!\n");
      exit(EXIT_FAILURE);
  }

  randomArr = calloc(NUMINDS*NUMLOCI,sizeof(double));
 
  /* ====== OpenCL initialized ====== */

  printf ("\n\n--------------------------------------\n\n");
  printf ("Finished initialization; starting MCMC \n");
  printf ("%d iterations + %d burnin\n\n", NUMREPS, BURNIN);

  /*=====Main MCMC loop=======================================*/

  for (rep = 0; rep < (NUMREPS + BURNIN); rep++) {

    UpdateP (P,LogP, Epsilon, Fst, NumAlleles, Geno, Z, lambda, Individual);

    if (LINKAGE && rep >= ADMBURNIN) {
      UpdateQMetroRecombine (Geno, Q, Z, P, Alpha, rep,
                             Individual, Mapdistance, R, Phase,Phasemodel);
    } else {
      UpdateQ (Geno, PreGeno, Q, P, Z, Alpha, rep, Individual, UsePopProbs, Recessive, LocPrior);
    }

    if (LOCPRIOR && UPDATELOCPRIOR) {
      UpdateLocPrior(Q, LocPrior, Alpha, Individual);
    }
    
    if (RECESSIVEALLELES) {
      UpdateGeno (PreGeno, Geno, P, Z, Recessive, NumAlleles, Q);
    /*The Zs are not correct after UpdateGeno, until UpdateZ is run */
    }
    
    if (LINKAGE && rep > ADMBURNIN) {
      if (!INDIVIDUALR) {
        recomblikelihood = UpdateZandSingleR(Z,  Q, P, Geno,
                                             R, Mapdistance, rep, Phase, Z1,Phasemodel, rep+1 > BURNIN? sumIndLikes : NULL, indLikesNorm);
      } else {
        recomblikelihood = UpdateZandR(Z,  Q, P, Geno, R,
                                       Mapdistance, rep, Phase, Z1,Phasemodel, rep+1 > BURNIN ? sumIndLikes:NULL, indLikesNorm);
      }
    } else {
      FillArrayWithRandom(randomArr,NUMINDS*NUMLOCI);
      compareZCLandZ(clDict,Z,Q,P,Geno,randomArr);
      UpdateZCL (clDict,Z,  Q, P, Geno,randomArr);
      /*
      UpdateZ (Z,  Q, P, Geno,randomArr);
      */
      /*      printf("done updatez alpha[2]=%e\n", Alpha[2]); */
    }

    if (LOCPRIOR && NOADMIX==0) {
      UpdateAlphaLocPrior(Q, Alpha, LocPrior, Individual);
    } else if (INFERALPHA) {
      UpdateAlpha (Q, Alpha, Individual, rep);
    }
    
    if (INFERLAMBDA) {
      if  (POPSPECIFICLAMBDA) {
        UpdatePopLambda(LogP,lambda,NumAlleles);
      } else {
        UpdateLambda (LogP,Epsilon,lambda, NumAlleles);
      }
    }


    if (FREQSCORR) {
      UpdateEpsilon(P,LogP,Epsilon,Fst,NumAlleles,lambda[0]);
      UpdateFst (Epsilon, Fst, LogP, NumAlleles);
    }

    /*====book-keeping stuff======================*/
    if (rep + 1 > BURNIN) {
      DataCollection (Geno, PreGeno, Q, QSum, Z, Z1,  P, PSum,
                      Fst, FstSum, NumAlleles,
                      AncestDist, Alpha, sumAlpha, sumR, varR, &like,
                      &sumlikes, &sumsqlikes, R, Epsilon,SumEpsilon,recomblikelihood,
                      lambda, sumlambda, Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes, indLikesNorm, rep);
    }
    
    if ((savefreq) && ((rep + 1) > BURNIN) && (((rep + 1 - BURNIN) % savefreq) == 0)
        && ((rep + 1) != NUMREPS + BURNIN)) {
      OutPutResults (Geno, rep + 1, savefreq, Individual, PSum, QSum,
                      FstSum, AncestDist, UsePopProbs, sumlikes,
                     sumsqlikes, sumAlpha, sumR, varR,
                     NumAlleles, Translation, 0, Markername, R,
                     SumEpsilon,
                     lambda,sumlambda,sumLocPrior, LocPriorLen,
                     sumIndLikes, indLikesNorm, argc,argv);
    }


    if (PRINTLIKES) {
      PrintLike (like, rep, Geno, PreGeno, Q, P,recomblikelihood);
    }
    
    if (((rep + 1) % UPDATEFREQ) == 0) {
      PrintUpdate (rep + 1, Geno, PreGeno, Alpha, Fst, P, Q, like,
                   sumlikes, sumsqlikes, NumAlleles, R, lambda,Individual,
                   recomblikelihood, Recessive, LocPrior, LocPriorLen);
    }
  }

  /*====final book-keeping====================================*/
  if ((rep % UPDATEFREQ) != 0) {
    PrintUpdate (rep, Geno, PreGeno, Alpha, Fst, P, Q, like, sumlikes,
                 sumsqlikes, NumAlleles,R, lambda, Individual,recomblikelihood,
                 Recessive, LocPrior, LocPriorLen);
  }

  OutPutResults (Geno, rep, savefreq, Individual, PSum, QSum,
                  FstSum, AncestDist, UsePopProbs,
                 sumlikes, sumsqlikes,
                 sumAlpha, sumR, varR, NumAlleles, Translation, 1,
                 Markername, R, SumEpsilon,
                 lambda,sumlambda,sumLocPrior, LocPriorLen,
                 sumIndLikes, indLikesNorm,
                 argc,argv);

  /*=====Closing everything down==============================*/
  FreeAll(Mapdistance, Phase, Phasemodel, lambda, sumlambda, Markername, Geno, PreGeno, Recessive,
            Individual, Translation, NumAlleles, Z, Z1, Q, P, LogP, R, sumR, varR, Epsilon, SumEpsilon,
            Fst, FstSum, NumLociPop, PSum, QSum,  AncestDist, UsePopProbs, LocPrior,
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm, clDict);
  free(randomArr);
  return (0);
}

/*==========================================

  Notes:


  ============================================*/
