#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"



/*------------------------------------------*/
/*Melissa added 7/12/07*/
void UpdateLocPriorNoAdmix(double *Q, double *LocPrior,
                           struct IND *Individual) {
  int i, j;
  double **dcount, newpp, diff;
  double logdiff, e1, e2, g1, g2, *eta;
  int loc, pop1, pop2;

  /* get counts in each location, cluster */
  dcount = malloc(NUMLOCATIONS*sizeof(double*));
  for (i=0; i<NUMLOCATIONS; i++) {
    dcount[i] = malloc(MAXPOPS*sizeof(double));
    for (j=0; j<MAXPOPS; j++) {
      dcount[i][j] = 0.0;
    }
  }
  for (i=0; i<NUMINDS; i++) {
    for (j=0; j<MAXPOPS; j++) {
      dcount[Individual[i].myloc][j] += (double)Q[QPos(i, j)];
    }
  }

  /* first update r (LocPrior[0]) */
  eta = &LocPrior[1];
  newpp = RandomReal(LocPrior[0]-LOCPRIORSTEP, LocPrior[0]+LOCPRIORSTEP);
  if (newpp > 0.0 && newpp < MAXLOCPRIOR) {
    logdiff = mylgamma(newpp) - mylgamma(LocPrior[0]);
    for (i=0; i<MAXPOPS; i++) {
      logdiff += mylgamma(LocPrior[0]*eta[i]) - mylgamma(newpp*eta[i]);
    }

    logdiff *= NUMLOCATIONS;
    for (i=0; i<NUMLOCATIONS; i++) {
      for (j=0; j<MAXPOPS; j++) {
        logdiff += (newpp-LocPrior[0])*eta[j]*log(LocPrior[LocPriorPos(i,j)]);
      }
    }

    if (logdiff >= 0.0 || RandomReal(0,1) < exp(logdiff)) {
      LocPrior[0] = newpp;
    }
  }

  /* now update eta */
  pop1 = RandomInteger(0, MAXPOPS-1);
  pop2 = RandomInteger(0, MAXPOPS-2);
  if (pop2>=pop1) {
    pop2++;
  }
  
  diff = RandomReal(0, ALPHAPROPSD);
  e1 = eta[pop1]+diff;
  e2 = eta[pop2]-diff;
  if (e1 < 1.0 && e2 >0.0) {
    /* don't need to consider prior on alpha if using Dirichlet(1,1,..,1) */
    logdiff = NUMLOCATIONS*(mylgamma(LocPrior[0]*eta[pop1]) + mylgamma(LocPrior[0]*eta[pop2])
                            -mylgamma(LocPrior[0]*e1)
                            -mylgamma(LocPrior[0]*e2));
    for (i=0; i<NUMLOCATIONS; i++) {
      logdiff += LocPrior[0]*((e1-eta[pop1])*log(LocPrior[LocPriorPos(i, pop1)]) +
                              (e2-eta[pop2])*log(LocPrior[LocPriorPos(i, pop2)]));
    }

    if (logdiff >= 0.0 || RandomReal(0, 1) < exp(logdiff)) {
      eta[pop1] = e1;
      eta[pop2] = e2;
    }
  }

  /* now update gamma */
  for (loc=0; loc<NUMLOCATIONS; loc++) {
    pop1 = RandomInteger(0, MAXPOPS-1);
    pop2 = RandomInteger(0, MAXPOPS-2);
    if (pop2 >= pop1) {
      pop2++;
    }
    diff = RandomReal(0, ALPHAPROPSD);
    g1 = LocPrior[LocPriorPos(loc, pop1)]+diff;
    g2 = LocPrior[LocPriorPos(loc, pop2)]-diff;
    if (g1 < 1.0 && g2 > 0.0) {
      logdiff =
          (LocPrior[0]*eta[pop1]-1.0)*(log(g1)-log(LocPrior[LocPriorPos(loc, pop1)])) +
          (LocPrior[0]*eta[pop2]-1.0)*(log(g2)-log(LocPrior[LocPriorPos(loc, pop2)])) +
          dcount[loc][pop1]*(log(g1)-log(LocPrior[LocPriorPos(loc, pop1)])) +
          dcount[loc][pop2]*(log(g2)-log(LocPrior[LocPriorPos(loc, pop2)]));
      if (logdiff >= 0.0 || RandomReal(0,1)  < exp(logdiff)) {
        LocPrior[LocPriorPos(loc, pop1)] = g1;
        LocPrior[LocPriorPos(loc, pop2)] = g2;
      }
    }
  }
  for (i=0; i<NUMLOCATIONS; i++) {
    free(dcount[i]);
  }
  free(dcount);
}

void UpdateLocPriorAdmix(double *Q, double *LocPrior, double *Alpha, struct IND *Individual) {
  double diff;
  double currPP, newPP, PPdiff, globalpha, localpha, alpharat;
  int k, loc;

  diff = 0.0;
  currPP = LocPrior[0];
  newPP = RandomReal(currPP - LOCPRIORSTEP, currPP + LOCPRIORSTEP);
  if (newPP > 0.0 && newPP < MAXLOCPRIOR) {
    PPdiff = newPP - currPP;
    diff = 0.0;
    for (k=0; k<MAXPOPS; k++) {
      globalpha = Alpha[k];
      for (loc=0; loc<NUMLOCATIONS; loc++) {
        localpha = Alpha[AlphaPos(loc, k)];
        alpharat = localpha/globalpha;
        diff += globalpha*PPdiff*log(localpha) -
                PPdiff*localpha -
                mylgamma(globalpha*newPP) +
                mylgamma(globalpha*currPP) +
                globalpha*newPP*log(newPP) -
                globalpha*currPP*log(currPP);
      }
    }
    if (diff > 0.0 || RandomReal(0, 1) < exp(diff)) {
      LocPrior[0] = newPP;
    }
  }
}

void UpdateLocPrior(double *Q, double *LocPrior, double *Alpha, struct IND *Individual) {
  if (NOADMIX) {
    UpdateLocPriorNoAdmix(Q, LocPrior, Individual);
  } else {
    UpdateLocPriorAdmix(Q, LocPrior, Alpha, Individual);
  }
}
