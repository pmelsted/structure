#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "datain.h"
#include "mymath.h"


/*---------------------------------------*/
double
Forward (int *Z, double *IndividualQ, double *P, int *Geno, double Rec, int ind,
         double *RTransitProb, double *Mapdistance, double *Phase,int *Phasemodel)
{
  long pop, pop2,  loc, line;
  double loglikelihood, temp, problinked, tempP00, tempP01,
      tempP10, tempP11,asum;
  double *sum1,*sum2;
  double sqrtunder = sqrt (UNDERFLO);
  sum1=calloc(MAXPOPS,sizeof(double));
  sum2=calloc(MAXPOPS,sizeof(double));
  if (sum1==NULL || sum2==NULL) {
    printf("WARNING: unable to allocate array space in Forward\n");
    Kill();
  }
  
  loglikelihood = 0.0;
  if (PHASED) {
    for (line = 0; line < LINES; line++) {
      /* set RTransitProb for the first locus, for each population */
      for (pop = 0; pop < MAXPOPS; pop++) {
        RTransitProb[RTransitProbPos (0, line, pop)] = IndividualQ[pop];
        if (Geno[GenPos (ind, line, 0)] != MISSING) {
          RTransitProb[RTransitProbPos (0, line, pop)] *= P[PPos (0, pop, Geno[GenPos (ind, line, 0)])];
        }
      }
      /* rest of the loci */
      for (loc = 1; loc < NUMLOCI; loc++) {
        if (Mapdistance[loc] < 0) {
          problinked=0;  /*this is the code for unlinked loci*/
        } else {
          problinked = exp (-Rec * Mapdistance[loc]);
        }
        
        /*these steps are time critical. Could save further time by storing values of exp(-Rec*mapdistance) */
        asum=0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          asum+=RTransitProb[RTransitProbPos(loc-1,line,pop)];
        }
        
        for (pop = 0; pop < MAXPOPS; pop++) {
          if (Geno[GenPos (ind, line, loc)] == MISSING) {
            RTransitProb[RTransitProbPos (loc, line, pop)] =
                (1.0-problinked)*IndividualQ[pop]*asum + problinked* RTransitProb[RTransitProbPos (loc - 1, line, pop)];
        /*necessary to incorporate recombination even if the data is missing */
          } else {
            RTransitProb[RTransitProbPos (loc, line, pop)] =
                ((1.0-problinked)*IndividualQ[pop]*asum + problinked * RTransitProb[RTransitProbPos (loc - 1, line, pop)])
                * P[PPos (loc, pop, Geno[GenPos (ind, line, loc)])];
          }
        }
        
        /*prevent underflows */
        temp = 0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          temp += RTransitProb[RTransitProbPos (loc, line, pop)];
        }
        if (temp < sqrtunder) {
          loglikelihood += log (sqrtunder);
          for (pop = 0; pop < MAXPOPS; pop++) {
            RTransitProb[RTransitProbPos (loc, line, pop)] /= sqrtunder;
          }
        }
      }

      temp = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++) {
        temp += RTransitProb[RTransitProbPos (NUMLOCI - 1, line, pop)];
      }
      loglikelihood += log (temp);
    }                           /* do the next LINE */
  } else {

    /* set RTransitProb for the first locus, for each population */
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
        if (Geno[GenPos (ind, 0, 0)] != MISSING) {
          tempP00 = P[PPos (0, pop, Geno[GenPos (ind, 0, 0)])];
          tempP01 = P[PPos (0, pop2, Geno[GenPos (ind, 0, 0)])];
        } else {
          tempP00 = 1.0;
          tempP01 = 1.0;
        }
        
        if (Geno[GenPos (ind, 1, 0)] != MISSING) {
          tempP10 = P[PPos (0, pop, Geno[GenPos (ind, 1, 0)])];
          tempP11 = P[PPos (0, pop2, Geno[GenPos (ind, 1, 0)])];
        } else {
          tempP10 = 1.0;
          tempP11 = 1.0;
        }
        
        if (Phasemodel[ind]==1) {
          RTransitProb[DiploidRTransitProbPos (0, pop, pop2)] =
              IndividualQ[pop] * IndividualQ[pop2] *
              (Phase[PhasePos (ind, 0)] * tempP00 * tempP11 + (1.0 - Phase[PhasePos (ind, 0)]) * tempP01 * tempP10);
        } else {
          RTransitProb[DiploidRTransitProbPos (0, pop, pop2)] = IndividualQ[pop] * IndividualQ[pop2] * tempP00 * tempP11;
        }
      }
    }

    /* rest of the loci */
    for (loc = 1; loc < NUMLOCI; loc++) {
      if (Mapdistance[loc] < 0) {
        problinked=0;  /*this is the code for unlinked loci*/
      } else {
        problinked = exp (-Rec * Mapdistance[loc]);
      }

      /*these steps are time critical. Could save further time by storing values of exp(-Rec*mapdistance) */
      /* first digit of pop indicates whether it refers to  locus loc (0) or locus loc-1 (loc)
         if MARKOVPHASE=0, second indicates whether it refers to maternal (0) or paternal (1)
         otherwise second indicates whether it is first (0) or second(1) locus  */
      for (pop=0;pop<MAXPOPS;pop++) {
        sum1[pop]=0.0;
        sum2[pop]=0.0;
      }
      
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
          sum1[pop]+=RTransitProb[DiploidRTransitProbPos(loc-1,pop,pop2)];
          sum2[pop2]+=RTransitProb[DiploidRTransitProbPos(loc-1,pop,pop2)];
        }
      }
      
      asum=0.0;
      for (pop=0;pop<MAXPOPS;pop++) {
        asum+=sum1[pop];
      }
      
      for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
        for (pop = 0; pop < MAXPOPS; pop++) {

          /* first digit of tempP indicates whether it refers to first (0) or second (1) listed allele copy.
             if Phasemodel=1
             Second digit indicates whether it refers to maternal (0) or paternal (1) strands
             otherwise
             second digit indicates whether it refers to first (0) or second (1) strands */

          if (Geno[GenPos (ind, 0, loc)] != MISSING) {
            tempP00 = P[PPos (loc, pop, Geno[GenPos (ind, 0, loc)])];
            tempP01 = P[PPos (loc, pop2, Geno[GenPos (ind, 0, loc)])];
          } else {
            tempP00 = 1.0;
            tempP01 = 1.0;
          }
          
          if (Geno[GenPos (ind, 1, loc)] != MISSING) {
            tempP10 = P[PPos (loc, pop, Geno[GenPos (ind, 1, loc)])];
            tempP11 = P[PPos (loc, pop2, Geno[GenPos (ind, 1, loc)])];
          } else {
            tempP10 = 1.0;
            tempP11 = 1.0;
          }

          /*note that for markov model (Phasemodel==0), phase information starts at locus 1  */
          if (Phasemodel[ind]==1) {
            RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)] = (
                problinked*problinked*RTransitProb[DiploidRTransitProbPos (loc - 1, pop, pop2)]
                + (1.0-problinked)*(1.0-problinked)*IndividualQ[pop]*IndividualQ[pop2]*asum
                +problinked*(1.0-problinked)*(IndividualQ[pop2]*sum1[pop]+IndividualQ[pop]*sum2[pop2])
                                                                      )* (Phase[PhasePos (ind, loc)] * tempP00 * tempP11 + (1.0 - Phase[PhasePos (ind, loc)]) * tempP01 * tempP10);
          } else {
            RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)] =  tempP00* tempP11*(
                problinked*problinked*
                (Phase[PhasePos(ind,loc)]* RTransitProb[DiploidRTransitProbPos(loc-1,pop,pop2)]+
                 (1.0-Phase[PhasePos(ind,loc)])*RTransitProb[DiploidRTransitProbPos(loc-1,pop2,pop)])
                +(1.0-problinked)*(1.0-problinked)*
                IndividualQ[pop]*IndividualQ[pop2]*asum
                +problinked*(1.0-problinked)*
                (Phase[PhasePos (ind, loc)]*(IndividualQ[pop2]*sum1[pop]+IndividualQ[pop]*sum2[pop2])
                 + (1.0-Phase[PhasePos (ind, loc)])*(IndividualQ[pop]*sum1[pop2]+IndividualQ[pop2]*sum2[pop]))
                                                                                        );
          }
        }
      }
      
      /*prevent underflows */
      temp = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++) {
        for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
          temp += RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)];
        }
      }

      if (temp < sqrtunder) {
        loglikelihood += log (sqrtunder);
        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)] /= sqrtunder;
          }
        }
      }
    }

    /* return the log of the sum of the RTransitProbs from the last locus, this time not summed over lines */
    temp = 0.0;
    for (pop = 0; pop < MAXPOPS; pop++) {
      for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
        temp += RTransitProb[DiploidRTransitProbPos (NUMLOCI - 1, pop, pop2)];
      }
    }
    loglikelihood += log (temp);
  }
  
  free(sum1);
  free(sum2);
  return loglikelihood;
}



/*-----------------------------------------------------------*/
void
Backward (int *Z, double *SiteBySiteSum, double *IndividualQ, double Rec, int ind,
          double *Mapdistance, double *RTransitProb, int rep, int *Z1,
          double *Phase, double *P, int *Geno,int *Phasemodel)
{
  int loc, pop, line, pop2, answer;
  double sum, sum2, problinked;
  double *Cutoffs, *Cutoffs2;
  double *SquareCutoffs;
  double temp00,temp01,temp10,temp11;
  /*  double temp; */

  /*added 1 in next two lines because the minimum allowable size is 2 (see last
    bit of this function where these start to refer to chromosome strands).*/
  Cutoffs = calloc (MAXPOPS+1, sizeof (double));
  Cutoffs2 = calloc (MAXPOPS+1, sizeof (double));
  SquareCutoffs = calloc (MAXPOPS * MAXPOPS, sizeof (double));
  if (Cutoffs == NULL || SquareCutoffs == NULL || Cutoffs2 == NULL) {
    printf ("WARNING: unable to allocate array space in Backwards\n");
    Kill ();
  }
  
  if (PHASED) {
    for (line = 0; line < LINES; line++)
    {
      /*NUMLOCI-1th locus */
      sum = 0.0;
      sum2 = 0.0;
      for (pop = 0; pop < MAXPOPS; pop++)
      {
        Cutoffs[pop] = RTransitProb[RTransitProbPos (NUMLOCI - 1, line, pop)];
        sum += Cutoffs[pop];
        Cutoffs2[pop] = P[PPos (NUMLOCI-1, pop, Geno[GenPos (ind, line, NUMLOCI-1)])];
        sum2 += Cutoffs2[pop];
      }

      Z[ZPos (ind, line, NUMLOCI - 1)] = PickAnOption (MAXPOPS, sum, Cutoffs);

      if (rep + 1 > BURNIN && SITEBYSITE) {
        for (pop = 0; pop < MAXPOPS; pop++)
          if (POSTERIOR)
            SiteBySiteSum[SiteBySiteSumPos (ind, line, NUMLOCI - 1, pop)] += Cutoffs[pop] / sum;
          else
            SiteBySiteSum[SiteBySiteSumPos (ind, line, NUMLOCI - 1, pop)] += Cutoffs2[pop] / sum2;
      }
      for (loc = NUMLOCI - 2; loc > -1; loc = loc - 1)
      {
        sum = 0.0;
        sum2 = 0.0;
        if (Mapdistance[loc+1] < 0) problinked=0;  /*this is the code for unlinked loci--JP*/
        else problinked = exp (-Rec * Mapdistance[loc + 1]);
        /*Same time saving approach as in Forward */
        /* here temp has bitten the dust */
        for (pop = 0; pop < MAXPOPS; pop++)
        {
          if (pop == Z[ZPos (ind, line, loc + 1)])
            Cutoffs[pop] = RTransitProb[RTransitProbPos (loc, line, pop)] * (problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, line, loc + 1)]]);
          else
            Cutoffs[pop] = RTransitProb[RTransitProbPos (loc, line, pop)] * (1.0 - problinked) * IndividualQ[Z[ZPos (ind, line, loc + 1)]];
          sum += Cutoffs[pop];
          Cutoffs2[pop] = P[PPos (loc, pop, Geno[GenPos (ind, line, loc)])];
          sum2 += Cutoffs2[pop];
        }
        if (rep + 1 > BURNIN && SITEBYSITE) {
          for (pop = 0; pop < MAXPOPS; pop++)
            if (POSTERIOR)
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] += Cutoffs[pop] / sum;
            else
              SiteBySiteSum[SiteBySiteSumPos (ind, line, loc, pop)] += Cutoffs2[pop] / sum2;
        }
        Z[ZPos (ind, line, loc)] = PickAnOption (MAXPOPS, sum, Cutoffs);
      }
    }
  } else {
    if (Phasemodel[ind]==1) {
      /*have treated first loci and others together to avoid excessive repetition */
      for (loc = NUMLOCI - 1; loc > -1; loc = loc - 1) {

        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            SquareCutoffs[SquarePos (pop, pop2)] = RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)];
          }
        }

        if (loc < NUMLOCI - 1) {
          if (Mapdistance[loc+1] < 0) {
            problinked=0;  /*this is the code for unlinked loci*/
          } else { problinked = exp (-Rec * Mapdistance[loc + 1]);
          }
          
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              if (pop == Z1[ZPos (ind, 0, loc + 1)]) {
                SquareCutoffs[SquarePos (pop, pop2)] *= problinked + (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 0, loc + 1)]];
              } else {
                SquareCutoffs[SquarePos (pop, pop2)] *= (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 0, loc + 1)]];
              }

              if (pop2 == Z1[ZPos (ind, 1, loc + 1)]) {
                SquareCutoffs[SquarePos (pop, pop2)] *= problinked + (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 1, loc + 1)]];
              } else {
                SquareCutoffs[SquarePos (pop, pop2)] *= (1.0 - problinked) * IndividualQ[Z1[ZPos (ind, 1, loc + 1)]];
              }
            }
          }
        }

        sum = 0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            sum += SquareCutoffs[SquarePos (pop, pop2)];
          }
        }

        if (rep + 1 > BURNIN && SITEBYSITE) {
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              SiteBySiteSum[DiploidSiteBySiteSumPos (ind, pop2, loc, pop)] += SquareCutoffs[SquarePos (pop, pop2)] / sum;
            }
          }
        }

        answer = PickAnOption (MAXPOPS * MAXPOPS, sum, SquareCutoffs);
        Z1[ZPos (ind, 0, loc)] = answer / MAXPOPS;
        Z1[ZPos (ind, 1, loc)] = answer - MAXPOPS * (int) (answer / MAXPOPS);
      }
      /* we have determined the populations for maternal and paternal strands Z1 */
      /*Now we work out the populations for the first and second loci Z */
      /*Note that meaning of Zpos changes from (ind,pop,loc) to (ind,line,loc). */
      for (loc = 0; loc < NUMLOCI; loc++) {
        Cutoffs[0] = (Phase[PhasePos(ind,loc)])
            * P[PPos (loc, Z1[ZPos (ind, 1, loc)], Geno[GenPos (ind, 1, loc)])]
            * P[PPos (loc, Z1[ZPos (ind, 0, loc)], Geno[GenPos (ind, 0, loc)])];
        Cutoffs[1] = (1.0 - Phase[PhasePos(ind,loc)])
            * P[PPos (loc, Z1[ZPos (ind, 0, loc)], Geno[GenPos (ind, 1, loc)])]
            * P[PPos (loc, Z1[ZPos (ind, 1, loc)], Geno[GenPos (ind, 0, loc)])];
        sum = Cutoffs[0] + Cutoffs[1];
        answer = PickAnOption (2, sum, Cutoffs);
        if (answer == 0) {
          Z[ZPos (ind, 0, loc)] = Z1[ZPos (ind, 0, loc)];
          Z[ZPos (ind, 1, loc)] = Z1[ZPos (ind, 1, loc)];
        } else {
          Z[ZPos (ind, 0, loc)] = Z1[ZPos (ind, 1, loc)];
          Z[ZPos (ind, 1, loc)] = Z1[ZPos (ind, 0, loc)];
        }
      }
    } else {
      /*have treated first loci and others together to avoid excessive repetition */
      for (loc = NUMLOCI - 1; loc > -1; loc = loc - 1) {

        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            SquareCutoffs[SquarePos (pop, pop2)] = RTransitProb[DiploidRTransitProbPos (loc, pop, pop2)];
          }
        }

        if (loc < NUMLOCI - 1) {
          if (Mapdistance[loc+1] < 0) {
            problinked=0;  /*this is the code for unlinked loci*/
          } else {
            problinked = exp (-Rec * Mapdistance[loc + 1]);
          }
          
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              if (pop == Z[ZPos (ind, 0, loc + 1)]) {
                temp00= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              } else {
                temp00= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              }

              if (pop2 == Z[ZPos (ind, 1, loc + 1)]) {
                temp11= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              } else {
                temp11= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              }

              if (pop2 == Z[ZPos (ind, 0, loc + 1)]) {
                temp01= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              } else {
                temp01= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 0, loc + 1)]];
              }

              if (pop == Z[ZPos (ind, 1, loc + 1)]) {
                temp10= problinked + (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              } else {
                temp10= (1.0 - problinked) * IndividualQ[Z[ZPos (ind, 1, loc + 1)]];
              }
              SquareCutoffs[SquarePos(pop,pop2)]*=temp00*temp11*Phase[PhasePos(ind,loc+1)]+temp10*temp01*(1.0-Phase[PhasePos(ind,loc+1)]);
            }
          }
        }

        sum = 0.0;
        for (pop = 0; pop < MAXPOPS; pop++) {
          for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
            sum += SquareCutoffs[SquarePos (pop, pop2)];
          }
        }

        if (rep + 1 > BURNIN && SITEBYSITE) {
          for (pop = 0; pop < MAXPOPS; pop++) {
            for (pop2 = 0; pop2 < MAXPOPS; pop2++) {
              SiteBySiteSum[ DiploidSiteBySiteSumPos (ind, pop2, loc, pop)] += SquareCutoffs[SquarePos (pop, pop2)] / sum;
            }
          }
        }

        answer = PickAnOption (MAXPOPS * MAXPOPS, sum, SquareCutoffs);
        Z[ZPos (ind, 0, loc)] = answer / MAXPOPS;
        Z[ZPos (ind, 1, loc)] = answer - MAXPOPS * (int) (answer / MAXPOPS);
      }
    }
  }

  free (Cutoffs);
  free (Cutoffs2);
  free (SquareCutoffs);
}
