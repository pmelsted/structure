#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran.h"
#include "mymath.h"
#include "structure.h"
#include "output.h"


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


