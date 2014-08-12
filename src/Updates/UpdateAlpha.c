#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"
#include "../Kernels.h"



/*-----------------------------------------*/
float AlphaPriorDiff (float newalpha, float oldalpha)
{
    /*returns log diff in priors for the alpha, assuming a gamma prior on alpha
      See notes 7/29/99 */
    return ((ALPHAPRIORA - 1) * log (newalpha / oldalpha) +
            (oldalpha - newalpha) / ALPHAPRIORB);
}

/*-----------------------------------------*/
float LogProbQDiff (float *Q, float oldalpha, float newalpha,
                     struct IND *Individual,float sumalphas,
                     int pop,int numredpops)
{
    /*return log prob of q given alpha [for single alpha in all populations].
      See notes 5/13/99 */
    float sum;
    int ind, redpop;
    /*this is the number of individuals without pop. info */
    int numinds = 0;
    float logterm = 0.0;
    int multiple = numredpops - pop;

    sum = 0.0;

    for (ind = 0; ind < NUMINDS; ind++) {
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            /* ie don't use individuals for whom prior pop info is used */
            numinds++;
            for (redpop = pop; redpop < numredpops; redpop++) {
                logterm += log(Q[QPos (ind, redpop)]);
            }
        }
    }
    sum -= (oldalpha - 1.0) * logterm;
    sum += (newalpha - 1.0) * logterm;

    sum -= (mylgamma (sumalphas) - multiple * mylgamma ( oldalpha)) * numinds;

    if (POPALPHAS){
        sumalphas += newalpha - oldalpha;
    } else {
        sumalphas = MAXPOPS*newalpha;
    }

    sum += (mylgamma (sumalphas) - multiple * mylgamma ( newalpha)) * numinds;

    return sum;
}

/*-----------------------------------------*/
float LogProbQTerm (float *Q, struct IND *Individual, int pop,int numredpops)
{
    /*return log prob of q given alpha [for single alpha in all populations].
      See notes 5/13/99 */
    int ind, redpop;
    /*this is the number of individuals without pop. info */
    float logterm = 0.0;


    for (ind = 0; ind < NUMINDS; ind++) {
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            /* ie don't use individuals for whom prior pop info is used */
            for (redpop = pop; redpop < numredpops; redpop++) {
                logterm += log(Q[QPos (ind, redpop)]);
            }
        }
    }
    return logterm;
}


/*-----------------------------------------*/
float LogProbQonepop (float *Q, float alpha, float alphasum,
                       struct IND *Individual,int pop)
{
    /*
     * return log prob of q given alpha--for one element of q ONLY.  This version is for
     * updates where there is one alpha for each population.  Everything cancels out of
     * the M-H ratio except one term of the gamma function, top and bottom, and the
     * relevant product in the q's
     */

    float sum;
    int ind;
    /*this is the number of individuals without pop. info */
    int numinds = 0;
    float logterm;

    sum = 0.0;
    logterm = 0.0;
    for (ind = 0; ind < NUMINDS; ind++) {
        if (!((USEPOPINFO) && (Individual[ind].PopFlag))) {
            /* ie don't use individuals for whom prior pop info is used */

            numinds++;
            /*being more careful with underflow caused by very small values of Q */
            /*0-values lead to underflow */
            logterm += log(Q[QPos (ind, pop)]);
        }
    }

    sum += (alpha - 1.0) * logterm;
    sum += (mylgamma (alphasum) - mylgamma (alpha)) * numinds;
    return (sum);
}

/* returns log Pr(Q|LocPrior) for a subset of individuals at a location */
float LogProbQ_LocPrior_loc(float *Q, float *Alpha, struct IND *Individual,
                             int loc)
{
    float sumalpha=0.0, sumgammaalpha=0.0, like=0.0;
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
    like += (mylgamma(sumalpha) - sumgammaalpha)*(float)numind;

    return like;
}

float LogProbQ (float *Q, float onealpha, struct IND *Individual)
{
  /*return log prob of q given alpha [for single alpha in all populations].
    See notes 5/13/99 */
  float sum;
  float runningtotal;
  int ind, pop;
  int numinds = 0;              /*this is the number of individuals without pop. info */
  float sqrtunder = sqrt (UNDERFLO);

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


void UpdateAlpha (float *Q, float *Alpha, struct IND *Individual, int rep)
{
  /* Produce new *Alpha using metropolis step.  There are two cases
     here: either there is the same alpha for all populations, or we do a
     separate Metropolis update for the alpha for each population.*/

  float newalpha;
  /*  float logoldprob;
   *  float lognewprob; */
  float unifrv;
  float threshold;
  float logprobdiff = 0;
  float sumalphas;
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


/* updates Alpha under LocPrior model */
void UpdateAlphaLocPrior(float *Q, float *Alpha, float *LocPrior,
                         struct IND *Individual)
{
    float diff, newalpha, oldalpha, lprobQ, globalpha, new_lprobQ;
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
            diff += (newalpha-oldalpha)*LocPrior[0]*log(Alpha[AlphaPos(loc,
                    pop)]) - mylgamma(newalpha*LocPrior[0]) + mylgamma(oldalpha*LocPrior[0]) +
                    (newalpha-oldalpha)*LocPrior[0]*log(LocPrior[0]);
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
            diff = (globalpha*LocPrior[0]-1.0)*log(newalpha/oldalpha) - LocPrior[0]*
                   (newalpha-oldalpha) + new_lprobQ - lprobQ;
            if (diff >= 0.0 || RandomReal(0,1) < exp(diff)) {
                lprobQ = new_lprobQ;
            } else {
                Alpha[pos+pop] = oldalpha;
            }
        }
    }
}

void UpdateAlphaCL (CLDict *clDict,float *Q, float *Alpha, struct IND *Individual, int rep, int POPFLAGINDS)
{
    /*
     * Produce new *Alpha using metropolis step.  There are two cases
     * here: either there is the same alpha for all populations, or we do a
     * separate Metropolis update for the alpha for each population.
     */

    /*float newalpha;*/
    /*float logprobdiff = 0;*/
    /*int pop, i;*/
    int numalphas;
    int numredpops;
    /*int redpop;*/
    size_t global[2];
    /*float alphasum;*/
    /*float oldalpha;*/
    /*float logterm;*/
    /*float sumalphas;*/
    /*float sum;*/
    /*int multiple;*/

    if (!((NOADMIX) && ((rep >= ADMBURNIN) || (rep > BURNIN)))) {
        /*don't update alpha in these cases*/
        if (POPALPHAS) {
            numredpops = 1;
            numalphas = MAXPOPS;
        } else {
            numredpops = MAXPOPS;
            numalphas = 1;
        }

        setKernelArg(clDict,PopNormals,ALPHACL,0);
        setKernelArgExplicit(clDict,PopNormals,sizeof(float),&ALPHAPROPSD,3);
        setKernelArgExplicit(clDict,PopNormals,sizeof(int),&numalphas,4);
        runTask(clDict,PopNormals,"PopNormals Alpha");

        /* The smaller the group, the less danger of underflow is there */
        /* global[0] = pow(2,(int) (log(NUMINDS)/log(2))); */
        global[0] = NUMINDS;
        global[1] = numalphas;
        runKernel(clDict,UpdateAlphaKernel,2,global,"Update Alpha kernel");
        /*
        alphasum = 0.0;
        for (i=0; i<MAXPOPS; i++)  {
            alphasum += Alpha[i];
        }

        for (pop = 0; pop < numalphas; pop++) {
            if (POPALPHAS){ numredpops = pop+1;}

            oldalpha = Alpha[pop];
            newalpha = RNormal (oldalpha, ALPHAPROPSD);
            [>reject immed. if out of range<]
            if ((newalpha > 0) && ((newalpha < ALPHAMAX) || (!(UNIFPRIORALPHA)) ) ) {
                if (!(UNIFPRIORALPHA)) {
                    logprobdiff = AlphaPriorDiff (newalpha, Alpha[pop]);
                }
                multiple = numredpops - pop;
                logterm = LogProbQTerm(Q,Individual,pop,numredpops);

                sum = 0.0;
                sum -= (oldalpha - 1.0) * logterm;
                sum += (newalpha - 1.0) * logterm;

                sumalphas = alphasum;
                sum -= (mylgamma (sumalphas) - multiple * mylgamma ( oldalpha)) * POPFLAGINDS;

                if (POPALPHAS){
                    sumalphas += newalpha - oldalpha;
                } else {
                    sumalphas = MAXPOPS*newalpha;
                }

                sum += (mylgamma (sumalphas) - multiple * mylgamma ( newalpha)) * POPFLAGINDS;
                logprobdiff += sum;
                [>c ompute probabilities <]
                if (rnd() < exp (logprobdiff)) {
                    for (redpop = pop; redpop < numredpops; redpop++) {
                        Alpha[redpop] = newalpha;
                    }
                }
            }
        }
        writeBuffer(clDict,Alpha,sizeof(float) * MAXPOPS,ALPHACL,"Alpha");*/
    }
}
