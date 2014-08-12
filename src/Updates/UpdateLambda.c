#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"

void UpdatePopLambda (float *P, float *lambda, int *NumAlleles)
/*updates a lambda for each population*/
{
    float new;
    float sum;
    float sumlogp;
    int loc, pop, allele;

    for (pop=0; pop<MAXPOPS; pop++) {
        new = RNormal (lambda[pop], LAMBDAPROPSD); /*proposal*/

        if ((new > 0.0) && (new < LAMBDAMAX)) {
            sum = 0.0;
            for (loc=0; loc < NUMLOCI; loc++) {   /*compute log of likelihood ratio*/
                if (NumAlleles[loc] > 1) {
                    /*norm constants*/
                    sum +=  mylgamma((float) NumAlleles[loc]*new);
                    sum -=  mylgamma((float) NumAlleles[loc]*lambda[pop]);
                    sum +=  (float) NumAlleles[loc] * mylgamma(lambda[pop]);
                    sum -=  (float) NumAlleles[loc] * mylgamma(new);

                    /*printf("%d %1.3f ----- ",loc,sum);*/

                    sumlogp = 0.0;
                    for (allele=0; allele<NumAlleles[loc]; allele++) {
                        sumlogp += log(P[PPos(loc,pop,allele)]);
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
void UpdateLambda (float *P,float *Epsilon, float *lambda,
                   int *NumAlleles)
/*
 * updates single value of lambda.  If FREQSCORR is turned on, this is based on the
 * ancestral frequencies (Epsilon); otherwise it is based on ALL
 * the population frequencies, P.  Uniform prior for lambda assumed
 */
{
    float new;
    float sum;
    float sumlogp;
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
                sum += (float) stoppop * mylgamma((float) NumAlleles[loc]*new);
                sum -= (float) stoppop * mylgamma((float) NumAlleles[loc]*lambda[0]);
                sum += (float) stoppop * (float) NumAlleles[loc] * mylgamma(lambda[0]);
                sum -= (float) stoppop * (float) NumAlleles[loc] * mylgamma(new);
                /*printf("%d %1.3f ----- ",loc,sum);*/

                sumlogp = 0.0;
                for (pop=0; pop<stoppop; pop++) {
                    for (allele=0; allele<NumAlleles[loc]; allele++) {
                        if (FREQSCORR) {
                            sumlogp += log(Epsilon[EpsPos(loc,allele)]);
                        } else {
                            sumlogp += log(P[PPos(loc,pop,allele)]);
                        }
                    }
                }
                sum += (new - lambda[0])*sumlogp;
                /*printf("%1.3f\n",sum);*/
            }
        }

        if (rnd() < exp(sum)) {
            for (pop=0; pop<MAXPOPS; pop++) {
                lambda[pop]=new;
            }
        }
    }
}
