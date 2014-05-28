#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran.h"
#include "mymath.h"
#include "datain.h"
#include "output.h"

#include "Kernels.h"


void InitializeZ (int *Geno, struct IND *Individual, int *Z);

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
                     double *sumIndLikes, double *indlike_norm, CLDict *clDict)

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

  InitCLDict(clDict);
}
