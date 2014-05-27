/*=======================================================

  STRUCTURE.C

  Program for inferring population structure using multilocus
  genotype data.

  Code written by Daniel Falush, Melissa Hubisz, and Jonathan Pritchard

  See additional details in README file.

  =========================================================*/
#define VERSION "2.3.4 (Jul 2012)"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "mymath.h"
#include "ran.h"
#include "params.h"
#include "datain.h"
#include "output.h"
#include "UpdateQ.h"
#include "UpdateZ.h"
#include "UpdateP.h"
#include "UpdateLocPrior.h"

void InitializeZ (int *Geno, struct IND *Individual, int *Z);

/*========================================================
  ==========================================================*/
void
Welcome (FILE * file)
{
  fprintf (file, "\n\n");
  fprintf (file, "----------------------------------------------------\n");
  fprintf (file, "STRUCTURE by Pritchard, Stephens and Donnelly (2000)\n");
  fprintf (file, "     and Falush, Stephens and Pritchard (2003)\n");
  fprintf (file, "       Code by Pritchard, Falush and Hubisz\n");
  fprintf (file, "             Version %s\n", VERSION);
  fprintf (file, "----------------------------------------------------\n");

  fprintf (file, "\n\n");
  fflush(file);
}
/*-------------------------------------------------------*/
void Kill ()                            /*exit program */
{
  printf ("\nExiting the program due to error(s) listed above.\n\n");
  exit (1);
}
/*---------------------------------------*/
void CheckParamCombinations()
{
  if ((LINES>2) && (USEPOPINFO==1))
  {
    printf(" USEPOPINFO does not work for ploidy >2\n");
    Kill();
  }
  if ((PHASEINFO) && (LINES!=2))
  {
    printf("phase info is only applicable to diploid data!! \n");
    Kill();
  }
  if (LINES ==2 && PHASEINFO==1 && PHASED ==0 && MARKOVPHASE==-1 && LINKAGE)
  {
    printf("You need to specify a phase model using the parameter MARKOVPHASE!! \n");
    Kill();
  }
  if (LINKAGE && !MAPDISTANCES)
  {
    printf("Map distance information is required in order to run the linkage model. \n");
    Kill();
  }

  if ((LINKAGE) && (!PHASED) && (LINES!=2))
  {
    printf("unphased data only permitted for diploids!! \n");
    Kill();
  }

  if ((LINKAGE) && (NOADMIX))
  {
    printf("please choose the LINKAGE option or the NOADMIX option, but not both \n");
    Kill();
  }

  if ((INFERLAMBDA) && (FREQSCORR))
  {
    printf("Warning!, choosing both INFERLAMBDA and FREQSCORR parameters may leave the computer with insufficient information to estimate either parameter accurately \n");
  }



  if (((NOADMIX) || (LINKAGE)) && ADMBURNIN >= BURNIN)
  {
    printf("The ADMBURNIN should finish before the BURNIN!! \n");
    Kill();
  }

  if ((PFROMPOPFLAGONLY) && (!(POPFLAG)))
  {
    printf("PFROMPOPFLAGONLY can only be turned on when the data file contains POPFLAG data\n");
    Kill();
  }
  if (LINES>2 && RECESSIVEALLELES && MISSING==NOTAMBIGUOUS)
  {
    printf("The code for missingdata (MISSING) should be set differently to the code (NOTAMBIGUOUS) for loci whose genotypes are known");
    Kill();
  }
  if (LOCPRIOR && LINKAGE) {
    printf("LOCPRIOR model is not compatible with linkage model\n");
    Kill();
  }
  if (LOCPRIOR && USEPOPINFO) {
    printf("LOCPRIOR model is not compatible with USEPOPINFO option\n");
    Kill();
  }
  /*  if (RANDOMIZE && SEED!=-1) {
      printf("Warning: Seed from parameter file will not be used as RANDOMIZE is set to 1.  SEED in output file will be random seed drawn from time\n");
      }*/  /* modified by JKP */

  if (RANDOMIZE)
    printf("Note: RANDOMIZE is set to 1. The random number generator will be initialized using the system clock, ignoring any specified value of SEED.\n");

}
/*---------------------------------------*/
/*void FreeAll (int *Geno, double *Mapdistance, char *Markername,
              struct IND *Individual, int *Translation,
              int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP, double *Epsilon,
              double *Fst, int *NumLociPop, double *PSum, double *QSum,
              double *FstSum, int *AncestDist, double *UsePopProbs, double *R,
              double *sumR, double *varR, double *LocPrior, double *sumLocPrior)
*/
void FreeAll(double *Mapdistance, double *Phase, int *Phasemodel, double *lambda, double *sumlambda,
             char *Markername, int *Geno, int* PreGeno, int* Recessive, struct IND *Individual,
             int *Translation, int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP,
             double *R, double *sumR, double *varR, double *Epsilon, double *SumEpsilon, double *Fst,
             double *FstSum, int *NumLociPop, double *PSum, double *QSum,
             int *AncestDist, double *UsePopProbs, double *LocPrior, double *sumLocPrior,
             double *Alpha, double *sumAlpha, double *sumIndLikes, double *indLikesNorm)
{
  /** these variables are calloc'd in main and freed in the same order */

  free (Mapdistance);
  
  free(Phase);
  if (LINES==2 && PHASED == 0) {
    free(Phasemodel);
  }
  


  free(lambda);
  free(sumlambda);

  
  free (Markername);
  free (Geno);

  if (RECESSIVEALLELES) {
    free(PreGeno);
    free(Recessive);
  }


  free (Individual);
  free (Translation);
  free (NumAlleles);
  
  free (Z);
  free (Z1);
  
  free (Q);
  free (P);
  
  free (LogP);
  free (R);
  free (sumR);
  free (varR);

  free (Epsilon);
  
  if (FREQSCORR) {
    free(SumEpsilon);
  }
  
  
  free (Fst);
  free (FstSum);
  free (NumLociPop);
  
  free (PSum);
  free (QSum);


  if (ANCESTDIST) {
    free (AncestDist);
  }

  if (USEPOPINFO)  {
    free (UsePopProbs);
  }

  if (LOCPRIOR) {
    free(LocPrior);
    free(sumLocPrior);
  }

  
  free(Alpha);
  free(sumAlpha);

  free(sumIndLikes);
  free(indLikesNorm);
  
}
/*---------------------------------------*/
void
PrintTranslation (int *Translation, int *NumAlleles)
{
  int loc, allele;
  printf ("Translation matrix:\n");
  for (loc = 0; loc < NUMLOCI; loc++)
  {
    for (allele = 0; allele < NumAlleles[loc]; allele++)
      printf ("%2d ", Translation[TransPos (loc, allele)]);
    printf ("\n");
  }
  printf ("\n");
  for (loc = 0; loc < NUMLOCI; loc++)
    printf ("%d ", NumAlleles[loc]);
  printf ("\n");


}

/*-------------------------------------------*/
void
InitFromGeogPop (int *Geno, struct IND *Individual, int *Z, int verbose)
{
  /*
   * initialize the population of origin of each allele. These are
   * started in their given populations (when there is population info).
   * It is assumed that the populations are numbered 1..K.  If the
   * population identifier is out of range, the allele is assigned to a
   * random population.  These is used to check that the MCMC is not
   * producing unexpected results because it couldn't find the mode.
   */

  int ind, line, loc;
  int poperrors = 0;
  int pop;

  if (!(POPDATA)) {
    InitializeZ (Geno, Individual, Z);
    if (verbose) {
      printf ("Starting from a random configuration because POP=0\n");
    }
  } else {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (line = 0; line < LINES; line++) {
        for (loc = 0; loc < NUMLOCI; loc++) {
          if (Geno[GenPos (ind, line, loc)] == MISSING) {
            Z[ZPos (ind, line, loc)] = UNASSIGNED;
          } else {
            pop = Individual[ind].Population;
            if ((pop > 0) && (pop <= MAXPOPS)) {
              Z[ZPos (ind, line, loc)] = pop - 1;
            } else {
              Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
              poperrors++;
            }
          }
        }
      }
    }
    if (verbose) {
      printf ("USING GIVEN POPULATION INFO TO SET INITIAL CONDITION\n");
    }
    if ((verbose) && (poperrors)) {
      printf ("WARNING: unable to initialize %d individuals to the predefined\n", poperrors);
      printf ("populations because their population names were not in the range\n");
      printf ("{1..%d}.  These individuals were initialized at random.\n",MAXPOPS);
    }
  }
}
/*---------------------------------------*/
void
InitializeZ (int *Geno, struct IND *Individual, int *Z)
{
  /* 
   * initialize the population of origin of each allele. I pick these at
   * random because this seems to produce better behaviour than starting
   * out with everybody in one pop, for instance.  I also set missing Data
   * to the unassigned pop from the start.  STARTATPOPINFO indicates that
   * individuals should be started at their given populations. It is
   * assumed that the populations are numbered 1..K.  If the population
   * identifier is out of range, the allele is assigned to a random
   * population.
   */

  int ind, line, loc;
  int allele;
  int pop;

  for (ind = 0; ind < NUMINDS; ind++) {
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        allele = Geno[GenPos (ind, line, loc)]; /*missing data */
        if (allele == MISSING) {
          Z[ZPos (ind, line, loc)] = UNASSIGNED;
        } else {  /*------data present-----------*/
          if ((STARTATPOPINFO) && (POPDATA)) {    /*use given pops as initial Z */
            pop = Individual[ind].Population;
            if ((pop > 0) && (pop <= MAXPOPS)) {
              Z[ZPos (ind, line, loc)] = pop - 1;
            } else {
              Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
            }
          } else {          /*initial Z random */
            Z[ZPos (ind, line, loc)] = RandomInteger (0, MAXPOPS - 1);
          }
        }
        /*printf("%d ",Z[ZPos(ind,line,loc)]); */
      }
      /*printf("\n"); */
    }
  }

  if ((STARTATPOPINFO) && (POPDATA)) {
    printf ("USING GIVEN POPULATION INFO TO SET INITIAL CONDITION\n");
  }
}
/*---------------------------------------*/
void
InitFreqPriors (double *Epsilon, double *Fst, int *Geno, int *NumAlleles)
{
  int ind, line, loc, allele,pop;
  int value;
  int *Count /*[MAXALLELES] */ ;        /*stores number of copies of each allele
                                          at present locus */
  int total;

  if (!(FREQSCORR)) {           /*allele frequencies uncorrelated across populations */
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Epsilon[EpsPos (loc, allele)] = LAMBDA;
      }
    }
    for (pop = 0; pop< MAXPOPS; pop++) {
      Fst[pop] = 0.5;
    }
  } else {                      /*correlated allele frequencies------------------- */
    Count = calloc (MAXALLELES, sizeof (int));
    if (Count == NULL) {
      printf ("Error in assigning memory, InitFreqPriors\n");
      Kill ();
    }

    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Count[allele] = 0;
      }
      total = 0;
      for (ind = 0; ind < NUMINDS; ind++) {
        for (line = 0; line < LINES; line++) {
          value = Geno[GenPos (ind, line, loc)];
          if (value != MISSING) {
            total++;
            if ((value < 0) || (value >= NumAlleles[loc])) {
              printf ("WARNING: expected allele value, InitFreqPriors: loc %d, allele %d\n", loc, value);
            } else {
              Count[value]++;
            }
          }
        }
      }

      /*Start Epsilon at (mean sample frequencies) */
      /* add lambda to all counts to ensure positive frequencies
       * for recessive model etc */
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Epsilon[EpsPos (loc, allele)] =
            (((double) LAMBDA +
              (double) Count[allele]) / ((double) NumAlleles[loc] *
                                         (double) LAMBDA +
                                         (double) total));
        /* ((double) Count[allele] / total); */
      }
      /*printf("\n"); */
    }


    for (pop= 0; pop < MAXPOPS; pop++) {
      Fst[pop] = FPRIORMEAN;    /*start Fst at the prior mean */
    }
    if (Count != NULL) {
      free (Count);
    }
  }                             /*end, correlated allele frequencies------------- */

}
/*---------------------------------------*/


void
CheckPopPriors (struct IND *Individual)
    /*
     * This function is called when USEPOPINFO==1 to check that the
     * prior information on populations is ok
     */
{
  int ind;
  int numnopop;

  if ((MIGRPRIOR < 0.0) || (MIGRPRIOR > 1.0)) {
    printf ("MIGRPRIOR (which is currently set to %1.3f) must be in the range [0.0, 1.0]\n", MIGRPRIOR);
    Kill ();
  }

  if (!(POPDATA)) {
    printf ("Can't apply USEPOPINFO because no POPDATA in input data file\n");
    Kill ();
  }

  if (!(POPFLAG)) {              /*if no popflag, then assume that everybody should be used */
    for (ind = 0; ind < NUMINDS; ind++) {
      Individual[ind].PopFlag = 1;
    }
  }

  /*Check that the given population is within range for all individuals, and if
    not, then turn popflag off for that individual. */
  for (ind = 0; ind < NUMINDS; ind++) {
    if (Individual[ind].PopFlag) {
      if ((Individual[ind].Population < 1) || (Individual[ind].Population > MAXPOPS)) {
        printf ("Warning: population prior for individual %d is %d, which is not\n",
                ind + 1, Individual[ind].Population);
        printf ("  in the range 1..%d.  Population prior for this individual will\n",
                MAXPOPS);
        printf ("  be ignored\n");
        Individual[ind].PopFlag = 0;
      }
    }
  }

  if ((INFERALPHA) && (!(NOADMIX))) {     /*check whether alpha is needed at all */
    numnopop = 0;
    for (ind = 0; ind < NUMINDS; ind++) {
      if (Individual[ind].PopFlag == 0) {
        numnopop++;
      }
    }
    if (numnopop == 0) {
      NOALPHA = 1;
      INFERALPHA = 0;
    }
  }
}


/*
 * GetNumLocations: Melissa added 7/12/07.  Sets the variable NUMLOCATIONS and also setse
 * all the individuals so that ind[i].myloc is in (0...NUMLOCATIONS).  ind[i].loc is unchanged
 * to whatever the input file indicates
 */
void GetNumLocations (struct IND *ind) {
  int maxloc=0, i, j, *freq, pos;
  for (i=0; i<NUMINDS; i++) {
    /* for now we're not dealing with unknown location */
    if (ind[i].Location < 0) {
      printf("individual %s has negative location!  locations should be >= 0\n", ind[i].Label);
      Kill();
    }
    if (ind[i].Location > maxloc) {
      maxloc = ind[i].Location;
    }
  }

  freq = malloc((maxloc+1)*sizeof(int));
  for (i=0; i<=maxloc; i++) {
    freq[i]=0;
  }
  for (i=0; i<NUMINDS; i++) {
    freq[ind[i].Location]++;
  }

  pos=0;
  for (i=0; i<=maxloc; i++) {
    if (freq[i]==0) {
      continue;
    }
    for (j=0; j<NUMINDS; j++) {
      if (ind[j].Location==i) {
        ind[j].myloc = pos;
      }
    }
    pos++;
  }
  free(freq);
  NUMLOCATIONS = pos;
}


/*---------------------------------------*/
void
InitializeSums (double *PSum, double *QSum,  double *FstSum,
                int *NumAlleles, int *AncestDist, double *UsePopProbs,double *SumEpsilon)
    /*initialize arrays which store sums of parameters */
{
  int loc, pop, allele, ind, box, gen, line;

  for (loc = 0; loc < NUMLOCI; loc++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        PSum[PPos (loc, pop, allele)] = 0.0;
      }
    }
  }

  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      QSum[QPos (ind, pop)] = 0.0;
    }
  }


  if (ANCESTDIST) {
    for (box = 0; box < NUMBOXES; box++) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (ind = 0; ind < NUMINDS; ind++) {
          AncestDist[AncestDistPos (ind, pop, box)] = 0;
        }
      }
    }
  }

  if (USEPOPINFO) {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (gen = 0; gen <= GENSBACK; gen++) {
          UsePopProbs[UsePPrPos (ind, pop, gen)] = 0.0;
        }
      }
    }
  }

  if (FREQSCORR) {
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        SumEpsilon[ EpsPos(loc,allele)] = 0.0;
      }
    }
    for (pop = 0; pop < MAXPOPS; pop++) {
      FstSum[pop] = 0.0;
    }
  }
}


/*-----------------------------------------------------*/
void
InitializeR (double *R, double *sumR, double *varR)
{
  int ind;
  for (ind = 0; ind < NUMINDS; ind++) {
    if (LOG10RSTART> LOG10RMIN && LOG10RSTART<LOG10RMAX) {
      R[ind]=exp(LOG10RSTART*2.302585092994046);
    } else {
      R[ind] = exp((LOG10RMIN+LOG10RMAX)*1.15129);
    }
  }
  sumR[ind] = 0.0;
  varR[ind] = 0.0;
}
/*---------------------------------------*/
void InitializeGeno (int *Geno, int *PreGeno)
{
  int ind, line, loc;
  for (ind = 0; ind < NUMINDS; ind++) {
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        Geno[GenPos (ind, line, loc)] = PreGeno[GenPos (ind, line, loc)];
      }
    }
  }
}

/*---------------------------------------*/
void Initialization (int *Geno, int *PreGeno,
                     struct IND *Individual, int *Translation,
                     int *NumAlleles, int *Z, int *Z1, double *Epsilon,
                     double *SumEpsilon,
                     double *Fst,double *PSum, double *Q, double *QSum,
                      double *FstSum,
                     int *AncestDist, double *UsePopProbs, double *Alpha,
                     double *sumAlpha, double *sumR, double *varR,
                     double *sumlikes, double *sumsqlikes,
                     int *savefreq, double *R, double *lambda, double *sumlambda,
                     double *Phase, int *Recessive, double *LocPrior,
                     double *sumLocPrior, int LocPriorLen,
                     double *sumIndLikes, double *indlike_norm)

    /*
     * This function is in charge of initializing the data arrays and other
     * parameters to appropriate values
     */
{
  int pop, ind, loc, i;

  NOALPHA = 0;   /*initialize the alphas*/
  for (pop=0; pop<MAXPOPS; pop++) {
    Alpha[pop] = ALPHA;
    sumAlpha[pop] = 0.0;
    lambda[pop]=LAMBDA;
    sumlambda[pop]=0.0;
  }

  if (LOCPRIOR) {
    if (NOADMIX==0) {
      for (loc=0; loc<NUMLOCATIONS; loc++) {
        for (pop=0; pop<MAXPOPS; pop++) {
          Alpha[AlphaPos(loc, pop)] = ALPHA;
          sumAlpha[AlphaPos(loc, pop)] = 0.0;
        }
      }
    }
    for (i=0; i<LocPriorLen; i++) {
      sumLocPrior[i]=0.0;
      LocPrior[i] = 1.0/(double)MAXPOPS;
    }
    LocPrior[0] = LOCPRIORINIT;
  }

  *sumlikes = 0.0;
  *sumsqlikes = 0.0;

  if (INTERMEDSAVE > 0) {
    *savefreq = (int) NUMREPS / (INTERMEDSAVE + 1);
    /*need to worry about saving one time too many due to integer truncation */
    if (((*savefreq) * (INTERMEDSAVE + 1)) < NUMREPS) {
      (*savefreq)++;
    }
  } else {
    *savefreq = 0;
  }

  if (LINKAGE) {
    InitializeR (R, sumR, varR);
  }

  if (RECESSIVEALLELES) {
    CountAlleles (PreGeno, NumAlleles, Translation, Recessive);
    InitializeGeno (Geno, PreGeno);
  } else {
    CountAlleles (Geno, NumAlleles, Translation, Recessive);  /*recode alleles to {0,..,1-k} */
  }

  if (STARTATPOPINFO) {
    InitFromGeogPop (Geno, Individual, Z, 1);    /*set Z to geog origin */
  } else {
    InitializeZ (Geno, Individual, Z);  /*set Z to random initial values */
  }

  InitFreqPriors (Epsilon, Fst, Geno, NumAlleles);      /*set priors on allele freqs */
  InitializeSums (PSum, QSum,  FstSum, NumAlleles, AncestDist, UsePopProbs,SumEpsilon);

  for (ind=0; ind<NUMINDS; ind++) {
    for (pop=0; pop<MAXPOPS; pop++) {
      Q[QPos(ind, pop)] = 1.0/MAXPOPS;
    }
  }

  if (USEPOPINFO) {
    CheckPopPriors (Individual);
    InitFromGeogPop (Geno, Individual, Z, 0);
  }

  /* bug -- corrected by Daniel April 2004
     if ((LINKAGE) && (!(PHASED)))
     for (ind = 0; ind < NUMINDS; ind++)
     for (loc = 0; loc < NUMLOCI; loc++)
     Phase[PhasePos (ind, loc)] = 0.5;
  */
  if ((LINKAGE) && (!(PHASED))&&(! (PHASEINFO))) {
    for (ind = 0; ind < NUMINDS; ind++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        Phase[PhasePos (ind, loc)] = 0.5;
      }
    }
  }

  for (ind=0; ind<NUMINDS; ind++) {
    sumIndLikes[ind] = indlike_norm[ind] = 0.0;
  }
}


/*-----------------------------------------*/
double LogProbQ (double *Q, double onealpha, struct IND *Individual)
{
  /*return log prob of q given alpha [for single alpha in all populations].
    See notes 5/13/99 */
  double sum;
  double runningtotal;
  int ind, pop;
  int numinds = 0;              /*this is the number of individuals without pop. info */
  double sqrtunder = sqrt (UNDERFLO);

  sum = 0.0;

  runningtotal = 1.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */
      numinds++;
      for (pop = 0; pop < MAXPOPS; pop++) {
        /*being more careful with underflow caused by very small values of Q */
        if (Q[QPos (ind, pop)] > sqrtunder) {             /*0-values lead to underflow */
          runningtotal *= Q[QPos (ind, pop)];
        } else {
          runningtotal *= sqrtunder;
        }

        if (runningtotal < sqrtunder)  {  /*this is to avoid having to take logs all the time */
          if (runningtotal == 0.0) {
            printf ("*");
          }
          sum += (onealpha - 1.0) * log (runningtotal);
          runningtotal = 1.0;
        }
      }
    }
  }

  sum += (onealpha - 1.0) * log (runningtotal);
  sum += (mylgamma (MAXPOPS * onealpha) - MAXPOPS * mylgamma (onealpha)) * numinds;

  /*printf("%1.2e ",sum);
    printf("%d ",numinds); */
  return (sum);
}


/*-----------------------------------------*/
double LogProbQonepop (double *Q, double popalpha, double sumalphas, struct IND *Individual,int pop)
{
  /*
   * return log prob of q given alpha--for one element of q ONLY.  This version is for
   * updates where there is one alpha for each population.  Everything cancels out of
   * the M-H ratio except one term of the gamma function, top and bottom, and the
   * relevant product in the q's
   */

  double sum;
  double runningtotal;
  int ind;
  int numinds = 0;              /*this is the number of individuals without pop. info */
  double sqrtunder = sqrt (UNDERFLO);
  sum = 0.0;
  runningtotal = 1.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
      /* ie don't use individuals for whom prior pop info is used */

      numinds++;
      /*being more careful with underflow caused by very small values of Q */
      /*0-values lead to underflow */

      if (Q[QPos (ind, pop)] > sqrtunder) {
        runningtotal *= Q[QPos (ind, pop)];
      } else {
        runningtotal *= sqrtunder;
      }

      if (runningtotal < sqrtunder) {    /*this is to avoid having to take logs all the time */
        if (runningtotal == 0.0) {
          printf ("*");
        }
        sum += (popalpha - 1.0) * log (runningtotal);
        runningtotal = 1.0;
      }
    }
  }

  sum += (popalpha - 1.0) * log (runningtotal);
  sum += (mylgamma (sumalphas) - mylgamma (popalpha)) * numinds;
  return (sum);
}


/*-----------------------------------------*/
double AlphaPriorDiff (double newalpha, double oldalpha)
{
  /*returns log diff in priors for the alpha, assuming a gamma prior on alpha
    See notes 7/29/99 */
  return ((ALPHAPRIORA - 1) * log (newalpha / oldalpha) + (oldalpha - newalpha) / ALPHAPRIORB);
}


/*-----------------------------------------*/
void UpdateAlpha (double *Q, double *Alpha, struct IND *Individual, int rep)
{
  /*
   * Produce new *Alpha using metropolis step.  There are two cases
   * here: either there is the same alpha for all populations, or we do a
   * separate Metropolis update for the alpha for each population.
   */

  double newalpha;
  /*  double logoldprob;
   *  double lognewprob; */
  double unifrv;
  double threshold;
  double logprobdiff = 0;
  double sumalphas;
  int pop, numalphas,i;

  if (!((NOADMIX) && ((rep >= ADMBURNIN) || (rep > BURNIN)))) {
    /*don't update alpha in these cases*/
    if (POPALPHAS) {
      numalphas = MAXPOPS;
    }
    else numalphas = 1;
    for (pop = 0; pop < numalphas; pop++) {
      newalpha = RNormal (Alpha[pop], ALPHAPROPSD); /*generate proposal alpha */

      /*reject immed. if out of range*/
      if ((newalpha > 0) && ((newalpha < ALPHAMAX) || (!(UNIFPRIORALPHA)) ) ) {
        if (!(UNIFPRIORALPHA)) {
          logprobdiff = AlphaPriorDiff (newalpha, Alpha[pop]);
        }
        /*compute probabilities */
        if (POPALPHAS)  { /*different alphas in each population*/
          sumalphas = 0.0;  /*need to send in sum of alphas*/
          for (i=0; i<MAXPOPS; i++)  {
            sumalphas += Alpha[i];
          }

          /*compute probabilities for M-H ratio*/
          logprobdiff -= LogProbQonepop (Q, Alpha[pop], sumalphas,Individual,pop);
          sumalphas += newalpha - Alpha[pop];
          logprobdiff += LogProbQonepop (Q, newalpha, sumalphas, Individual,pop);
        } else  {  /*same alpha for all populations*/
          logprobdiff += LogProbQ (Q, newalpha, Individual);
          logprobdiff -= LogProbQ (Q, Alpha[pop], Individual);
        }

        /*accept new alpha with min of 1 and exp(logprobdiff) */
        threshold = exp (logprobdiff);
        unifrv = rnd ();

        /*printf("%d %.3f %.3f %.4f     ",pop,Alpha[pop],newalpha,threshold);
          if (pop==MAXPOPS-1) printf("\n");*/

        if (unifrv < threshold) {
          Alpha[pop] = newalpha;

          if (!(POPALPHAS)) { /*if same alpha in all populations*/
            for (pop = 1; pop < MAXPOPS; pop++) {
              Alpha[pop] = newalpha;
            }
          }
        }
      }
    }
  }
}

/*--------------------------------------------*/

/* returns log(exp(a)+exp(b)) without causing overflow issues */
double logsumexp(double a, double b)
{
  if (a-b > 100) {
    return a;
  }
  if (b-a > 100) {
    return b;
  }
  return a + log(1.0 + exp(b-a));
}

/* returns log Pr(Q|LocPrior) for a subset of individuals at a location */
double LogProbQ_LocPrior_loc(double *Q, double *Alpha, struct IND *Individual, int loc)
{
  double sumalpha=0.0, sumgammaalpha=0.0, like=0.0;
  int ind, pop, numind=0;

  for (ind=0; ind<NUMINDS; ind++) {
    if (Individual[ind].myloc!=loc) {
      continue;
    }
    for (pop=0; pop<MAXPOPS; pop++) {
      like += (Alpha[pop]-1.0)*log(Q[QPos(ind, pop)]);
    }
    numind++;
  }

  for (pop=0; pop<MAXPOPS; pop++) {
    sumalpha += Alpha[pop];
    sumgammaalpha += mylgamma(Alpha[pop]);
  }
  like += (mylgamma(sumalpha) - sumgammaalpha)*(double)numind;

  return like;
}


/* updates Alpha under LocPrior model */
void UpdateAlphaLocPrior(double *Q, double *Alpha, double *LocPrior,
                         struct IND *Individual)
{
  double diff, newalpha, oldalpha, lprobQ, globalpha, new_lprobQ;
  int pop, loc, pos;

  /* first update global alpha */
  for (pop=0; pop < MAXPOPS; pop++) {
    oldalpha = Alpha[pop];
    newalpha = RNormal(oldalpha, ALPHAPROPSD);
    if (newalpha >= ALPHAMAX || newalpha <= 0.0) {
      continue;
    }
    diff = 0.0;
    for (loc=0; loc<NUMLOCATIONS; loc++) {
      diff += (newalpha-oldalpha)*LocPrior[0]*log(Alpha[AlphaPos(loc,pop)]) - mylgamma(newalpha*LocPrior[0]) + mylgamma(oldalpha*LocPrior[0]) + (newalpha-oldalpha)*LocPrior[0]*log(LocPrior[0]);
    }

    if (diff > 0.0 || RandomReal(0,1) < exp(diff)) {
      Alpha[pop] = newalpha;
    }
  }

  /* now update location-specific alphas */
  for (loc=0; loc<NUMLOCATIONS; loc++) {
    pos = AlphaPos(loc, 0);
    lprobQ = LogProbQ_LocPrior_loc(Q, &Alpha[pos], Individual, loc);
    for (pop=0; pop<MAXPOPS; pop++) {
      globalpha = Alpha[pop];
      oldalpha = Alpha[pos+pop];
      newalpha = RNormal(oldalpha, ALPHAPROPSD);

      if (newalpha <= 0.0) {
        continue;
      }
      Alpha[pos+pop] = newalpha;
      new_lprobQ = LogProbQ_LocPrior_loc(Q, &Alpha[pos], Individual, loc);
      diff = (globalpha*LocPrior[0]-1.0)*log(newalpha/oldalpha) - LocPrior[0]*(newalpha-oldalpha) + new_lprobQ - lprobQ;
      if (diff >= 0.0 || RandomReal(0,1) < exp(diff)) {
        lprobQ = new_lprobQ;
      } else {
        Alpha[pos+pop] = oldalpha;
      }
    }
  }
}


/*--------------------------------------------*/

void UpdatePopLambda (double *LogP, double *lambda, int *NumAlleles)
    /*updates a lambda for each population*/
{
  double new;
  double sum;
  double sumlogp;
  int loc, pop, allele;

  for (pop=0;pop<MAXPOPS;pop++) {
    new = RNormal (lambda[pop], LAMBDAPROPSD); /*proposal*/

    if ((new > 0.0) && (new < LAMBDAMAX)) {
      sum = 0.0;
      for (loc=0; loc < NUMLOCI; loc++) {   /*compute log of likelihood ratio*/
        if (NumAlleles[loc] > 1) {
          /*norm constants*/
          sum +=  mylgamma((double) NumAlleles[loc]*new);
          sum -=  mylgamma((double) NumAlleles[loc]*lambda[pop]);
          sum +=  (double) NumAlleles[loc] * mylgamma(lambda[pop]);
          sum -=  (double) NumAlleles[loc] * mylgamma(new);

          /*printf("%d %1.3f ----- ",loc,sum);*/

          sumlogp = 0.0;
          for (allele=0; allele<NumAlleles[loc]; allele++) {
            sumlogp += LogP[PPos(loc,pop,allele)];
          }
          sum += (new - lambda[pop])*sumlogp;
          /*printf("%1.3f\n",sum);*/
        }
      }

      if (rnd() < exp(sum)) {
        lambda[pop]=new;
      }
    }
  }
}


/*--------------------------------------------*/
void UpdateLambda (double *LogP,double *Epsilon, double *lambda, int *NumAlleles)
    /*
     * updates single value of lambda.  If FREQSCORR is turned on, this is based on the
     * ancestral frequencies (Epsilon); otherwise it is based on ALL
     * the population frequencies, P.  Uniform prior for lambda assumed
     */
{
  double new;
  double sum;
  double sumlogp;
  int loc, pop, allele,stoppop;

  new = RNormal (lambda[0], LAMBDAPROPSD); /*proposal*/

  if ((new > 0.0) && (new < LAMBDAMAX)) {
    if (FREQSCORR) {
      stoppop=1;
    } else {
      stoppop = MAXPOPS;
    }

    sum = 0.0;
    for (loc=0; loc < NUMLOCI; loc++) {  /*compute log of likelihood ratio*/
      if (NumAlleles[loc] > 1) {
        /*norm constants*/
        sum += (double) stoppop * mylgamma((double) NumAlleles[loc]*new);
        sum -= (double) stoppop * mylgamma((double) NumAlleles[loc]*lambda[0]);
        sum += (double) stoppop * (double) NumAlleles[loc] * mylgamma(lambda[0]);
        sum -= (double) stoppop * (double) NumAlleles[loc] * mylgamma(new);
        /*printf("%d %1.3f ----- ",loc,sum);*/

        sumlogp = 0.0;
        for (pop=0; pop<stoppop; pop++) {
          for (allele=0; allele<NumAlleles[loc]; allele++) {
            if (FREQSCORR) {
              sumlogp += log(Epsilon[EpsPos(loc,allele)]);
            } else {
              sumlogp += LogP[PPos(loc,pop,allele)];
            }
          }
        }
        sum += (new - lambda[0])*sumlogp;
        /*printf("%1.3f\n",sum);*/
      }
    }

    if (rnd() < exp(sum)) {
      for (pop=0;pop<MAXPOPS;pop++) {
        lambda[pop]=new;
      }
    }
  }
}


/*============================================*/
double
FPriorDiff (double newf, double oldf)
{
  /*returns log diff in priors for the correlation, f. See notes 5/14/99, and 7/15/99 */

  return ((FPRIORMEAN*FPRIORMEAN/(FPRIORSD*FPRIORSD) - 1) * log (newf / oldf) + (oldf - newf) *FPRIORMEAN/(FPRIORSD*FPRIORSD));

}


/*-----------------------------------------*/
double
FlikeFreqs (double f, double *Epsilon, double *LogP, int *NumAlleles, int pop)
{
  /*
   * returns the log probability of the allele frequencies (for a particular pop)
   * given the prior and the z (number of each allele in each population).
   * Here f is the value of Fst for that locus
   */

  /*
   * If numalleles=1 this seems to be ok
   * here passes Epsilon into mylgamma. Does this cause problems if epsilon very small?
   */
  int allele;
  double sum;
  int loc;
  double frac = (1.0-f)/f;

  sum = NUMLOCI*mylgamma(frac);
  for (loc=0; loc<NUMLOCI; loc++) {
    for (allele=0; allele < NumAlleles[loc]; allele++) {
      sum += frac*Epsilon[EpsPos (loc, allele)]*LogP[PPos(loc,pop,allele)];
      sum -= mylgamma( frac*Epsilon[EpsPos (loc, allele)]);
    }
    if (NumAlleles[loc]==0) {
      sum -=mylgamma(frac); /* should not be counting sites with all missing data */
    }
  }
  return sum;
}


/*-----------------------------------------*/
void
UpdateFst (double *Epsilon, double *Fst,
           double *LogP, int *NumAlleles)
    /*update the correlation factor, Fst, for each population*/
{

  double newf,oldf;
  double logprobdiff;
  double unifrv;
  double threshold;
  int pop1,pop2;
  int numpops1, numpops2;

  /*------Update f ()----See notebook, 5/14/99-----------*/

  /*There are two models: either there is a different F for each population,
    in which case we go through the entire loop K times; otherwise there
    is a single F, in which case we sum the likelihood ratio across populations.*/

  /*control the outer loop*/
  if (ONEFST) {
    numpops1 = 1;
  } else {
    numpops1 = MAXPOPS;
  }

  for (pop1 = 0; pop1 < numpops1; pop1++) {
    /*generate proposal f */
    oldf = Fst[pop1];
    newf = RNormal (oldf, FPRIORSD);

    /*reject if propopal < 0 or greater than 1 */
    if (newf > 0.0 && newf<1.0) {
      /*compute prior ratio */
      logprobdiff = FPriorDiff (newf, oldf);

      /*compute log likelihood diff */
      if (ONEFST) {
        numpops2 = MAXPOPS;
      } else {
        numpops2 = pop1+1;
      }
      for (pop2 = pop1; pop2 < numpops2; pop2++){
        logprobdiff += FlikeFreqs (newf, Epsilon, LogP, NumAlleles, pop2);
        logprobdiff -= FlikeFreqs (oldf, Epsilon, LogP, NumAlleles, pop2);
      }

      /*decide whether to accept, and then update*/

      if (logprobdiff >= 0.0) {   /*accept new f */
        for (pop2 = pop1; pop2 < numpops2; pop2++) {
          Fst[pop2] = newf;
        }
      } else {                 /*accept new parameter with prob p */
        threshold = exp (logprobdiff);
        unifrv = rnd ();
        if (unifrv < threshold) {
          for (pop2 = pop1; pop2 < numpops2; pop2++) {
            Fst[pop2] = newf;
          }
        }
      }
    }
  }
}






/*------------------------------------------*/
void IndependenceUpdateEpsilon(double *P,double *LogP, double *Epsilon, double *Fst,int *NumAlleles, double Lambda)
    /*this is the alternative update to the one below, proposed by Graham */
{
  int loc, pop, allele;
  /*  double difference; */
  double Sum;
  double frac;
  double *trialepsilon,*parameters;

  trialepsilon = calloc (MAXALLELES, sizeof (double));
  parameters=calloc(MAXALLELES,sizeof(double));

  if (trialepsilon == NULL || parameters == NULL) {
    printf ("warning: unable to allocate memory in UpdateEpsilon\n");
    Kill ();
  }

  for (loc = 0; loc < NUMLOCI; loc++) {
    if (NumAlleles[loc] > 1) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        parameters[allele] = Lambda;
        for (pop = 0; pop < MAXPOPS; pop++) {
          parameters[allele] +=
              (1.0 - Fst[pop]) * P[PPos (loc, pop, allele)] / Fst[pop];
        }
      }

      RDirichlet (parameters, NumAlleles[loc], trialepsilon);
      Sum = 0.0;
      /* compute Hastings (transition) ratio and prior ratio*/
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Sum +=
            (parameters[allele] - Lambda) * (log(Epsilon[EpsPos (loc, allele)]) -
                                             log(trialepsilon[allele]));
      }

      /*compute likelihood ratio*/

      for (pop = 0; pop < MAXPOPS; pop++) {
        frac = (1.0 - Fst[pop]) / Fst[pop];
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
          Sum +=
              mylgamma (frac * Epsilon[EpsPos (loc, allele)]);
          Sum -=
              mylgamma (frac * trialepsilon[allele]);
          Sum +=
              frac * (trialepsilon[allele] - Epsilon[EpsPos (loc, allele)])
              * LogP[PPos (loc,pop,allele)];
        }
      }

      if (rnd () < exp (Sum)) {
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
          Epsilon[EpsPos (loc, allele)] = trialepsilon[allele];
        }
      }
    }
  }
  free (trialepsilon);
  free (parameters);
}



/*------------------------------------------*/
void
UpdateEpsilon(double *P,double *LogP, double *Epsilon, double *Fst,int *NumAlleles, double lambda)
    /*
     * update the ancestral allele freq vector Epsilon.  This is done
     * by picking 2 alleles at each locus, and changing their frequencies.
     * See notes May 30; June 20, 2001
     */

{
  int loc,pop,allele1,allele2;
  double difference,invsqrtnuminds;
  double sum;
  double frac;
  /* double ratio; */

  /*here we choose between two different updates that we believe have different mixing
    properties, especially for small lambda. The independence update uses a
    Dirichlet prior independent of current epsilon while the update below uses a small normal jump */
  if (rnd()<0.5) {
    IndependenceUpdateEpsilon(P,LogP, Epsilon, Fst,NumAlleles, lambda);
  } else {
    /*this sets the range from which the proposal is drawn*/
    invsqrtnuminds=pow((double)NUMINDS,-0.5);

    for (loc=0;loc<NUMLOCI;loc++) {
      if (NumAlleles[loc]>1) {
        allele1=RandomInteger(0,NumAlleles[loc]-1);

        do {
          allele2=RandomInteger(0,NumAlleles[loc]-1);
        } while (allele1==allele2);

        difference=RandomReal(0,invsqrtnuminds);

        /*check that the proposals are in range*/
        if ((Epsilon[EpsPos(loc,allele1)]+difference<1.0) &&
            (Epsilon[EpsPos(loc,allele2)]-difference>0.0)) {

          sum=0.0;
          for (pop=0; pop<MAXPOPS; pop++) { /*compute likelihood ratio*/
            frac = (1.0-Fst[pop])/Fst[pop];

            sum += mylgamma(frac*Epsilon[EpsPos (loc, allele1)]);
            sum += mylgamma(frac*Epsilon[EpsPos (loc, allele2)]);
            sum -= mylgamma(frac*(Epsilon[EpsPos (loc, allele1)]+difference));
            sum -= mylgamma(frac*(Epsilon[EpsPos (loc, allele2)]-difference));

            sum += frac*difference*LogP[PPos (loc, pop, allele1)];
            sum -= frac*difference*LogP[PPos (loc, pop, allele2)];
          }

          if (lambda != 1.0) {              /*compute prior ratio*/
            /*TEMP: log added by JKP 6/30/03 as I think this was previously
              an error.  Now doing testing */
            sum += log(pow( (Epsilon[EpsPos (loc, allele1)] + difference)*
                            (Epsilon[EpsPos (loc, allele2)] - difference)/
                            (Epsilon[EpsPos (loc, allele1)])/
                            (Epsilon[EpsPos (loc, allele2)]), (double) lambda-1.0));
          }

          /*if (loc==3)
            {
            printf("%1.3f->%1.3f   %1.3f->%1.3f     ",Epsilon[EpsPos(loc,0)],
            Epsilon[EpsPos(loc,0)]
            +(allele1==0)*difference-(allele1==1)*difference,
            Epsilon[EpsPos(loc,1)],
            Epsilon[EpsPos(loc,1)]
            +(allele2==0)*difference-(allele2==1)*difference);
            printf("%1.3f %1.3f     MH=%1.5f\n",
            P[PPos (loc, 0, 0)],
            P[PPos (loc, 1, 0)],
            exp(sum));
            }*/


          if (rnd() < exp(sum)) {
            Epsilon[EpsPos(loc,allele1)]+=difference;
            Epsilon[EpsPos(loc,allele2)]-=difference;
          }
        }
      }
    }
  }
}



/*------------------------------------------*/
void UpdateGeno (int *PreGeno, int *Geno, double *P, int *Z,
                 int *Recessive, int *NumAlleles, double *Q)
    /* 
     * this function updates the imputed genotypes when the genotypes are
     * ambiguous due to recessive markers or inbreeding.
     */
{

  int ind, loc, allele, line, dom, toggle,rejectioncount,allelecount,notmissingcount;
  int *AlleleUsed=NULL, *AllelePresent=NULL;
  double *AlleleProbs, Sum;
  static int RejectionThreshold = 1000000;
  /*  int pop;
   *  double temp, Sum1, Sum2; */

  if (LINES == 2) {
    AlleleProbs = calloc (4, sizeof (double));
  } else {
    AlleleProbs = calloc (MAXALLELES, sizeof (double));
    AlleleUsed = calloc (MAXALLELES, sizeof (int));
    AllelePresent = calloc (MAXALLELES, sizeof (int));
  }

  if (LINES == 2) {  /* special update for diploids; general case by rejection sampling */
    for (ind = 0; ind < NUMINDS; ind++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        if (PreGeno[GenPos (ind, 0, loc)] != MISSING
            && PreGeno[GenPos (ind, 1, loc)] != MISSING) {
          if (PreGeno[GenPos (ind, 0, loc)] ==
              PreGeno[GenPos (ind, 1, loc)]) {
            for (dom = 0; dom < 4; dom++) {
              AlleleProbs[dom] = 0.0;
            }

            if (RECESSIVEALLELES
                && Recessive[loc] != MISSING    /* bug fixed 05072007 */
                && PreGeno[GenPos (ind, 0, loc)] != Recessive[loc]) {
              AlleleProbs[0] =
                  P[PPos(loc, Z[ZPos(ind, 0, loc)], Recessive[loc])] *
                  P[PPos(loc, Z[ZPos(ind, 1, loc)], PreGeno[GenPos(ind, 0, loc)])];
              AlleleProbs[1] =
                  P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                  P[PPos(loc, Z[ZPos (ind, 1, loc)], Recessive[loc])];
            }

            AlleleProbs[2] =
                P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                P[PPos(loc, Z[ZPos (ind, 1, loc)], PreGeno[GenPos (ind, 1, loc)])];

            Sum = AlleleProbs[0] + AlleleProbs[1] + AlleleProbs[2];
            dom = PickAnOption (3, Sum, AlleleProbs);

            if (dom == 0) {
              Geno[GenPos (ind, 0, loc)] = Recessive[loc];
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else if (dom == 1) {
              Geno[GenPos (ind, 1, loc)] = Recessive[loc];
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else if (dom == 2) {
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else {
              printf ("surely some mistake in UpdateGeno\n");
              Kill(); /* modify by William - kill is not a standard ANSI C function */
            }
          }
        }
      }
    }
  } else { /* general n-ploid case.  Select genotype by rejection sampling */
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (ind = 0; ind < NUMINDS; ind++) {
        if (Recessive[loc]==NOTAMBIGUOUS) {
          for (line=0;line<LINES;line++) {
            Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
          }
        } else {
          allelecount=0;
          notmissingcount=0;
          
          for (allele = 0; allele < NumAlleles[loc]; allele++) {
            AllelePresent[allele] = 0;
          }
          
          for (line = 0; line < LINES; line++) {
            if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
              notmissingcount+=1;
              if (AllelePresent[PreGeno[GenPos (ind, line, loc)]]==0) {
                AllelePresent[PreGeno[GenPos (ind, line, loc)]] = 1;
                allelecount+=1;
              }
            }
          }
          
          if (allelecount==notmissingcount) {  /* if number of alleles equal to number of slots then nothing to do */
            for (line=0;line<LINES;line++) {
              Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
            }
          } else {
            toggle = 1;
            rejectioncount=0;
            while (toggle && rejectioncount < RejectionThreshold) {
              rejectioncount+=1;
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                AlleleUsed[allele] = 0;
              }
              
              for (line = 0; line < LINES; line++) {
                if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
                  for (allele = 0; allele < NumAlleles[loc]; allele++) {
                    if (AllelePresent[allele] || allele == Recessive[loc]) {
                      AlleleProbs[allele] =
                          P[PPos (loc, Z[ZPos (ind, line, loc)], allele)];
                    } else {
                      AlleleProbs[allele] = 0.0;
                    }
                  }
                  Sum = 0.0;
                  for (allele = 0; allele < NumAlleles[loc]; allele++) {
                    Sum += AlleleProbs[allele];
                  }
                  dom = PickAnOption (NumAlleles[loc], Sum, AlleleProbs);
                  Geno[GenPos (ind, line, loc)] = dom;
                  AlleleUsed[dom] = 1;
                }
              }
              
              toggle = 0;
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (AlleleUsed[allele] == 0 && AllelePresent[allele] == 1) {
                  toggle = 1;
                }
              }
            }
            
            if (toggle==1) {
              /* rejection method failed, set a lower rejection threshold for future iterations */
              if (RejectionThreshold > 100) {
                RejectionThreshold /= 2;
                if (RejectionThreshold <100) {
                  RejectionThreshold = 100;
                }
              } 
              /*
                printf("sorry, STRUCTURE has tried to simulate alleles for individual %d at locus %d 1,000,000 times by a rejection method and failed", ind+1, loc+1);
                Kill();
              */
              line=0;

              /*  first fill in one copy of each present allele. */
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (AllelePresent[allele] == 1) {
                  Geno[GenPos(ind,line,loc)]= allele;
                  line++;
                }
              }

              /* now fill in remaining nonmissing slots. */
              for (line=allelecount; line<notmissingcount; line++) {
                for (allele = 0; allele < NumAlleles[loc]; allele++) {
                  if (AllelePresent[allele] || allele == Recessive[loc]) {
                    AlleleProbs[allele] = P[PPos (loc, Z[ZPos (ind, line, loc)], allele)];
                  } else {
                    AlleleProbs[allele] = 0.0;
                  }
                }

                Sum = 0.0;
                for (allele = 0; allele < NumAlleles[loc]; allele++) {
                  Sum += AlleleProbs[allele];
                }
                dom = PickAnOption (NumAlleles[loc], Sum, AlleleProbs);
                Geno[GenPos (ind, line, loc)] = dom;
                AlleleUsed[dom] = 1;
              }

              /* now fill in missing slots */
              for (line=notmissingcount;line<LINES;line++) {
                Geno[GenPos (ind, line, loc)]=MISSING;
              }
            }
          }
        }
      }
    }
  }
  free (AlleleProbs);
  if (LINES != 2) {
    free (AlleleUsed);
    free (AllelePresent);
  }
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
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm);
    Kill ();
  }
  /*=========done setting aside memory space=====================*/

  /*initialize variables and arrays */
  Initialization (Geno, PreGeno, Individual, Translation, NumAlleles, Z, Z1, Epsilon, SumEpsilon,
                  Fst, PSum, Q, QSum, FstSum, AncestDist, UsePopProbs, Alpha,
                  sumAlpha, sumR, varR, &sumlikes, &sumsqlikes, &savefreq, R, lambda,
                  sumlambda,Phase,Recessive, LocPrior, sumLocPrior, LocPriorLen, sumIndLikes, indLikesNorm);
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
      UpdateZ (Z,  Q, P, Geno,rep);
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
            sumLocPrior, Alpha, sumAlpha, sumIndLikes, indLikesNorm);
  return (0);
}

/*==========================================

  Notes:


  ============================================*/
