#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"
#include "../Kernels.h"

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
FlikeFreqs (double f, double *Epsilon, double *P, int *NumAlleles, int pop)
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
                sum += frac*Epsilon[EpsPos (loc, allele)]* log(P[PPos(loc,pop,allele)]);
                sum -= mylgamma( frac*Epsilon[EpsPos (loc, allele)]);
            }
        }
    }
    return sum;
}

double FlikeFreqsDiffMap (double newfrac,double oldfrac, double *Epsilon, double *P, int *NumAlleles, int loc,int pop){
    int allele;
    double eps,logp;
    double sum;

    if (NumAlleles[loc]==0) {
        return -(mylgamma(newfrac) - mylgamma(oldfrac)); /* should not be counting sites with all missing data */
    } else {
        sum = 0.0;
        for (allele=0; allele < NumAlleles[loc]; allele++) {
            eps = Epsilon[EpsPos (loc, allele)];
            logp = log(P[PPos(loc,pop,allele)]);
            sum += (newfrac-oldfrac)*eps*logp - (mylgamma( newfrac*eps) - mylgamma( oldfrac*eps));
        }
        return sum;
    }
}

/*-----------------------------------------*/
double
FlikeFreqsDiff (double newf,double oldf, double *Epsilon, double *P, int *NumAlleles, int pop)
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
    int loc;
    double newfrac = (1.0-newf)/newf;
    double oldfrac = (1.0-oldf)/oldf;
    double sum = 0.0;
    double init = NUMLOCI*(mylgamma(newfrac) - mylgamma(oldfrac));

    for (loc=0; loc<NUMLOCI; loc++) {
        sum += FlikeFreqsDiffMap(newfrac,oldfrac,Epsilon,P,NumAlleles,loc,pop);
    }
    return sum+init;
}



/*-----------------------------------------*/
void
UpdateFstCL (CLDict *clDict,double *Epsilon, double *Fst, double *P, int *NumAlleles)
/*update the correlation factor, Fst, for each population*/
{

    size_t global[2];
    int numfst;
    /*double *reduceresult;
    static int rep = 0;

    reduceresult = calloc(MAXGROUPS*NUMINDS*NUMLOCI,sizeof(double));
    writeBuffer(clDict,reduceresult,sizeof(double)*MAXGROUPS*NUMINDS*NUMLOCI,REDUCERESULTSCL,"result");*/

    /*------Update f ()----See notebook, 5/14/99-----------*/

    /*There are two models: either there is a different F for each population,
      in which case we go through the entire loop K times; otherwise there
      is a single F, in which case we sum the likelihood ratio across populations.*/

    /*control the outer loop*/
    numfst = MAXPOPS;
    if (ONEFST){
        numfst = 1;
    }

    global[0] = numfst;
    setKernelArg(clDict,PopNormals,FSTCL,0);
    setKernelArgExplicit(clDict,PopNormals,sizeof(double),&FPRIORSD,3);
    runKernel(clDict,PopNormals,1,global,"PopNormals Fst");

    global[0] = pow(2,(int) (log(NUMLOCI)/log(2)));
    /*global[0] = (128 < NUMLOCI ) ? 128 : NUMLOCI;*/
    global[1] = numfst;
    runKernel(clDict,UpdateFstKernel,2,global,"UpdateFst");
    /*if (rep % 100 == 0){
        readBuffer(clDict,reduceresult,sizeof(double)*MAXGROUPS*NUMINDS*NUMLOCI,REDUCERESULTSCL,"result");
        printf("%.2f %.2f \n",reduceresult[0],reduceresult[1]);
    }
    free(reduceresult);
    rep +=1; */
}
