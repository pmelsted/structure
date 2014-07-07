#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"

/*============================================*/
double
FPriorDiff (double newf, double oldf)
{
    /*returns log diff in priors for the correlation, f. See notes 5/14/99, and 7/15/99 */

    return ((FPRIORMEAN*FPRIORMEAN/(FPRIORSD*FPRIORSD) - 1) * log (newf / oldf) +
            (oldf - newf) *FPRIORMEAN/(FPRIORSD*FPRIORSD));

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
        if (NumAlleles[loc]==0) {
            sum -=mylgamma(frac); /* should not be counting sites with all missing data */
        } else {
            for (allele=0; allele < NumAlleles[loc]; allele++) {
                sum += frac*Epsilon[EpsPos (loc, allele)]*LogP[PPos(loc,pop,allele)];
                sum -= mylgamma( frac*Epsilon[EpsPos (loc, allele)]);
            }
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
            for (pop2 = pop1; pop2 < numpops2; pop2++) {
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
