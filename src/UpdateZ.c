#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran.h"
#include "mymath.h"
#include "structure.h"
#include "output.h"
#include "ForwardAndBackward.h"
/*-----------------------------------------*/
/*O*(NUMINDS*LINES*NUMLOCI*MAXPOPS)*/
void
UpdateZ (int *Z, double *SiteBySiteSum, double *Q, double *P, int *Geno,int rep)
    /*update Z: population origin of each allele */
{
  int ind, line, loc, pop;
  double *Cutoffs, *Cutoffs2 /*[MAXPOPS] */ ;
  /*Cutoffs contains unnormalized probabilities of
    an allele coming from each population */
  /*Cutoffs2 shows the same thing, independent of the individual's other alleles*/
  double sum=0.0, sum2=0.0;
  int allele;

  Cutoffs = calloc (MAXPOPS, sizeof (double));
  Cutoffs2 = calloc (MAXPOPS, sizeof (double));
  if (Cutoffs == NULL || Cutoffs2 == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZ\n");
    Kill ();
  }
  /*O*(NUMINDS*LINES*NUMLOCI*MAXPOPS)*/
  for (ind = 0; ind < NUMINDS; ind++) {  /*go through all alleles in sample */
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        allele = Geno[GenPos (ind, line, loc)];

        if (allele == MISSING) {   /*Missing Data */
          Z[ZPos (ind, line, loc)] = UNASSIGNED;
        } else {
          /*Data present */
          sum = 0.0;    /*compute prob of each allele being from each pop */
          sum2 = 0.0;
          for (pop = 0; pop < MAXPOPS; pop++) {
            Cutoffs[pop] = Q[QPos (ind, pop)] * P[PPos (loc, pop, allele)];
            sum += Cutoffs[pop];
          }

          Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs);
        }
        
      }
    }
  }

  free (Cutoffs);
  free (Cutoffs2);
}






/*----------------------------------------*/
double
UpdateZandSingleR (int *Z, double *SiteBySiteSum, double *Q, double *P, int *Geno,
                   double *R, double *Mapdistance, int rep, double *Phase,
                   int *Z1,int *Phasemodel, double *sumIndLikes,
                   double *indlike_norm)
    /* updates Z and R, assuming that the data is phased */
{
  long ind, pop;
  double *RTransitProb;
  double *IndividualQ;
  double trialR, logtrialR,currentloglikelihood, trialloglikelihood, indlike;
  /* long loc; */
  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  IndividualQ = calloc (MAXPOPS, sizeof (double));
  /*this form ensures compatibility with UpdateQMetroRecombine */

  if (RTransitProb == NULL || IndividualQ == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZandSingleR\n");
    Kill ();
  }
  
  currentloglikelihood = 0.0;
  trialloglikelihood = 0.0;
  logtrialR = RNormal(log(R[0])/2.30259,LOG10RPROPSD);
  if (logtrialR<LOG10RMIN) {
    logtrialR=2*LOG10RMIN-logtrialR;
  }

  if (logtrialR>LOG10RMAX) {
    logtrialR=2*LOG10RMAX-logtrialR;
  }
  
  trialR=exp(2.30259*logtrialR);
  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      IndividualQ[pop] = Q[QPos (ind, pop)];
    }
    indlike = Forward(Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                      Mapdistance, Phase,Phasemodel);
    currentloglikelihood += indlike;
    if (sumIndLikes!=NULL) {
      sumIndLikes[ind] += exp(indlike-indlike_norm[ind]);
    }

    Backward(Z, SiteBySiteSum, IndividualQ, R[ind], ind, Mapdistance,
             RTransitProb, rep, Z1, Phase, P, Geno,Phasemodel);

    trialloglikelihood += Forward (Z, IndividualQ, P, Geno, trialR, ind, RTransitProb,
                                   Mapdistance, Phase,Phasemodel);

  }
  /*printf("% .3E % 3E % .8E % .8E % .5E\n",trialR,R[0],trialloglikelihood,currentloglikelihood,exp (trialloglikelihood - currentloglikelihood)); */
  if (RandomReal (0.0, 1.0) < exp (trialloglikelihood - currentloglikelihood)) {  /*Accept */
    R[0] = trialR;
    /*currentloglikelihood=trialloglikelihood;  commented out by JKP--see email from Daniel 9 Dec 02*/
  }
  
  for (ind = 0; ind < NUMINDS; ind++) {
    R[ind] = R[0];
  }
  
  free (RTransitProb);
  free (IndividualQ);

  return currentloglikelihood;
}



/*----------------------------------------*/
double
UpdateZandR (int *Z, double *SiteBySiteSum, double *Q, double *P, int *Geno,
             double *R, double *Mapdistance, int rep, double *Phase, int *Z1,int *Phasemodel, double *sumindlike, double *indlike_norm)
    /* updates Z and R, assuming that the data is phased */
{
  long ind, pop;
  double *RTransitProb;
  double *IndividualQ;
  double trialR,logtrialR, currentloglikelihood, trialloglikelihood,sumlikelihood;
  /*  long loc; */


  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  IndividualQ = calloc (MAXPOPS, sizeof (double));
  /*this form ensures compatibility with UpdateQMetroRecombine */

  if (RTransitProb == NULL || IndividualQ == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZandR\n");
    Kill ();
  }
  
  sumlikelihood=0.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      IndividualQ[pop] = Q[QPos (ind, pop)];
    }

    currentloglikelihood = Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                                    Mapdistance, Phase,Phasemodel);
    Backward (Z, SiteBySiteSum, IndividualQ, R[ind], ind, Mapdistance, RTransitProb,
              rep, Z1, Phase, P, Geno,Phasemodel);
    
    logtrialR = RNormal(log(R[ind])/2.30259,LOG10RPROPSD);
    if (logtrialR<LOG10RMIN) {
      logtrialR=2*LOG10RMIN-logtrialR;
    }
    if (logtrialR>LOG10RMAX) {
      logtrialR=2*LOG10RMAX-logtrialR;
    }
    trialR=exp(2.30259*logtrialR);
    trialloglikelihood = Forward (Z, IndividualQ, P, Geno, trialR, ind, RTransitProb, Mapdistance, Phase,Phasemodel);
    /*printf("% .3E % 3E % .8E % .8E % .5E\n",trialR,R[ind],trialloglikelihood,currentloglikelihood,exp (trialloglikelihood - currentloglikelihood)); */
    if (RandomReal (0.0, 1.0) < exp (trialloglikelihood - currentloglikelihood)) {        /*Accept */
      R[ind] = trialR;
      sumlikelihood+=trialloglikelihood;
      if (sumindlike!=NULL) {
        sumindlike[ind] += exp(trialloglikelihood-indlike_norm[ind]);
      }
    } else {
      sumlikelihood+=currentloglikelihood;
      if (sumindlike!=NULL) {
        sumindlike[ind] += exp(currentloglikelihood-indlike_norm[ind]);
      }
    }
  }
  free (RTransitProb);
  free (IndividualQ);

  return sumlikelihood;
}
