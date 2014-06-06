#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"



/*===============================================*/
/*O(MAXPOPS*numalleles + NUMINDS*LINES) */
void
GetNumFromPop (int *NumAFromPop, int *Geno, int *Z, int loc,
	       int numalleles,struct IND *Individual)
{
  /*Fill in the number of each allele from each pop */
  int ind, line, pop, allele;
  /* int genpos; */
  int allelevalue;
  int popvalue;

  /* O(MAXPOPS*numalleles) */
  for (pop = 0; pop < MAXPOPS; pop++) {
    for (allele = 0; allele < numalleles; allele++) {
      NumAFromPop[NumAFromPopPos (pop, allele)] = 0;
    }
  }

  /* O(NUMINDS*LINES) */
  if (PFROMPOPFLAGONLY) {     /*this option uses only individuals with POPFLAG=1 to update P*/
    for (ind = 0; ind < NUMINDS; ind++) {
      if (Individual[ind].PopFlag == 1) {    /*individual must have popflag turned on*/
        for (line = 0; line < LINES; line++) {
          popvalue = Z[ZPos (ind, line, loc)];
          allelevalue = Geno[GenPos (ind, line, loc)];

          if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
            NumAFromPop[NumAFromPopPos (popvalue, allelevalue)]++;
          }
        }
      }
    }
  } else {       /*standard update--use everybody to update P */
    for (ind = 0; ind < NUMINDS; ind++) {
      for (line = 0; line < LINES; line++) {
        popvalue = Z[ZPos (ind, line, loc)];
        allelevalue = Geno[GenPos (ind, line, loc)];

        if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
          NumAFromPop[NumAFromPopPos (popvalue, allelevalue)]++;
        }
      }
    }
  }
}

/*------------------------------------------*/
/*
 * O(NUMLOCI*(MAXPOPS* (max_loc NumAlleles[loc]) + NUMINDS*LINES)) =>
 * O(NUMLOCI*(MAXPOPS*MAXALLELES + NUMINDS*LINES))
 */
void UpdateP (double *P, double *LogP, double *Epsilon, double *Fst,
              int *NumAlleles, int *Geno, int *Z, double *lambda, struct IND *Individual)
    /*Simulate new allele frequencies from Dirichlet distribution */
{
  int loc, pop, allele;
  double *Parameters;           /*[MAXALLS] **Parameters of posterior on P */
  int *NumAFromPop;             /*[MAXPOPS][MAXALLS] **number of each allele from each pop */

  Parameters = calloc(MAXALLELES, sizeof (double));
  NumAFromPop = calloc(MAXPOPS * MAXALLELES, sizeof (int));

  if ((Parameters == NULL) || (NumAFromPop == NULL)) {
    printf ("WARNING: unable to allocate array space in UpdateP\n");
    Kill ();
  }

  /* O(NUMLOCI*(MAXPOPS* (max_loc NumAlleles[loc]) + NUMINDS*LINES)) */
  for (loc = 0; loc < NUMLOCI; loc++) {
    /*count number of each allele from each pop */
    /*O(MAXPOPS*NumAlleles[loc] + NUMINDS*LINES) */
    GetNumFromPop (NumAFromPop, Geno, Z, loc, NumAlleles[loc], Individual);
    /* O(MAXPOPS*NumAlleles[loc])*/
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        if (FREQSCORR) {
          Parameters[allele] = Epsilon[EpsPos (loc, allele)]
              *(1.0- Fst[pop])/Fst[pop]
              + NumAFromPop[NumAFromPopPos (pop, allele)];
        } else {
          Parameters[allele] = lambda[pop]
              + NumAFromPop[NumAFromPopPos (pop, allele)];
        }
      }
      /*return a value of P simulated from the posterior Di(Parameters) */
      /*O(NumAlleles[loc]) */
      LogRDirichlet (Parameters, NumAlleles[loc],
		     P + PPos (loc, pop, 0),
		     LogP +PPos(loc,pop,0));

      /*need to worry about underflow in UpdateEpsilon due to
        allele frequencies being set to zero---hence previously used the
        following hack, however now pass LogP instead

        for (allele=0;allele<NumAlleles[loc];allele++) if
        (P[PPos(loc,pop,allele)]<1E-20) {
        P[PPos(loc,pop,allele)]=1E-20;

        for (pop=0; pop<MAXPOPS; pop++)
        {
        printf(" loc =%d pop= %d fst=%f ",loc,pop,Fst[pop]);
        for (allele=0;allele<NumAlleles[loc];allele++)
        printf (" Epsilon= %.5E P= %.5E Parameters=  %.5E Num= %d",
        Epsilon[EpsPos(loc,allele)],P[PPos(loc,pop,allele)],
        Parameters[allele],NumAFromPop[NumAFromPopPos (pop, allele)]);
        printf("\n");
        }
        } */
    }
  }

  free (Parameters);
  free (NumAFromPop);
}
