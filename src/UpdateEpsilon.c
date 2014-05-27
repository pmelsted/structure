#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran.h"
#include "mymath.h"
#include "structure.h"
#include "output.h"


/*------------------------------------------*/
void IndependenceUpdateEpsilon(double *P,double *LogP, double *Epsilon, double *Fst,int *NumAlleles, double Lambda)
    /*this is the alternative update to the one below, proposed by Graham */
{
  int loc, pop, allele;
  /*  double difference; */
  double Sum;
  double frac;
  double *trialepsilon,*parameters;

  trialepsilon = calloc (MAXALLELES, sizeof (double));
  parameters=calloc(MAXALLELES,sizeof(double));

  if (trialepsilon == NULL || parameters == NULL) {
    printf ("warning: unable to allocate memory in UpdateEpsilon\n");
    Kill ();
  }

  for (loc = 0; loc < NUMLOCI; loc++) {
    if (NumAlleles[loc] > 1) {
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        parameters[allele] = Lambda;
        for (pop = 0; pop < MAXPOPS; pop++) {
          parameters[allele] +=
              (1.0 - Fst[pop]) * P[PPos (loc, pop, allele)] / Fst[pop];
        }
      }

      RDirichlet (parameters, NumAlleles[loc], trialepsilon);
      Sum = 0.0;
      /* compute Hastings (transition) ratio and prior ratio*/
      for (allele = 0; allele < NumAlleles[loc]; allele++) {
        Sum +=
            (parameters[allele] - Lambda) * (log(Epsilon[EpsPos (loc, allele)]) -
                                             log(trialepsilon[allele]));
      }

      /*compute likelihood ratio*/

      for (pop = 0; pop < MAXPOPS; pop++) {
        frac = (1.0 - Fst[pop]) / Fst[pop];
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
          Sum +=
              mylgamma (frac * Epsilon[EpsPos (loc, allele)]);
          Sum -=
              mylgamma (frac * trialepsilon[allele]);
          Sum +=
              frac * (trialepsilon[allele] - Epsilon[EpsPos (loc, allele)])
              * LogP[PPos (loc,pop,allele)];
        }
      }

      if (rnd () < exp (Sum)) {
        for (allele = 0; allele < NumAlleles[loc]; allele++) {
          Epsilon[EpsPos (loc, allele)] = trialepsilon[allele];
        }
      }
    }
  }
  free (trialepsilon);
  free (parameters);
}

/*------------------------------------------*/
void
UpdateEpsilon(double *P,double *LogP, double *Epsilon, double *Fst,int *NumAlleles, double lambda)
    /*
     * update the ancestral allele freq vector Epsilon.  This is done
     * by picking 2 alleles at each locus, and changing their frequencies.
     * See notes May 30; June 20, 2001
     */

{
  int loc,pop,allele1,allele2;
  double difference,invsqrtnuminds;
  double sum;
  double frac;
  /* double ratio; */

  /*here we choose between two different updates that we believe have different mixing
    properties, especially for small lambda. The independence update uses a
    Dirichlet prior independent of current epsilon while the update below uses a small normal jump */
  if (rnd()<0.5) {
    IndependenceUpdateEpsilon(P,LogP, Epsilon, Fst,NumAlleles, lambda);
  } else {
    /*this sets the range from which the proposal is drawn*/
    invsqrtnuminds=pow((double)NUMINDS,-0.5);

    for (loc=0;loc<NUMLOCI;loc++) {
      if (NumAlleles[loc]>1) {
        allele1=RandomInteger(0,NumAlleles[loc]-1);

        do {
          allele2=RandomInteger(0,NumAlleles[loc]-1);
        } while (allele1==allele2);

        difference=RandomReal(0,invsqrtnuminds);

        /*check that the proposals are in range*/
        if ((Epsilon[EpsPos(loc,allele1)]+difference<1.0) &&
            (Epsilon[EpsPos(loc,allele2)]-difference>0.0)) {

          sum=0.0;
          for (pop=0; pop<MAXPOPS; pop++) { /*compute likelihood ratio*/
            frac = (1.0-Fst[pop])/Fst[pop];

            sum += mylgamma(frac*Epsilon[EpsPos (loc, allele1)]);
            sum += mylgamma(frac*Epsilon[EpsPos (loc, allele2)]);
            sum -= mylgamma(frac*(Epsilon[EpsPos (loc, allele1)]+difference));
            sum -= mylgamma(frac*(Epsilon[EpsPos (loc, allele2)]-difference));

            sum += frac*difference*LogP[PPos (loc, pop, allele1)];
            sum -= frac*difference*LogP[PPos (loc, pop, allele2)];
          }

          if (lambda != 1.0) {              /*compute prior ratio*/
            /*TEMP: log added by JKP 6/30/03 as I think this was previously
              an error.  Now doing testing */
            sum += log(pow( (Epsilon[EpsPos (loc, allele1)] + difference)*
                            (Epsilon[EpsPos (loc, allele2)] - difference)/
                            (Epsilon[EpsPos (loc, allele1)])/
                            (Epsilon[EpsPos (loc, allele2)]), (double) lambda-1.0));
          }

          /*if (loc==3)
            {
            printf("%1.3f->%1.3f   %1.3f->%1.3f     ",Epsilon[EpsPos(loc,0)],
            Epsilon[EpsPos(loc,0)]
            +(allele1==0)*difference-(allele1==1)*difference,
            Epsilon[EpsPos(loc,1)],
            Epsilon[EpsPos(loc,1)]
            +(allele2==0)*difference-(allele2==1)*difference);
            printf("%1.3f %1.3f     MH=%1.5f\n",
            P[PPos (loc, 0, 0)],
            P[PPos (loc, 1, 0)],
            exp(sum));
            }*/


          if (rnd() < exp(sum)) {
            Epsilon[EpsPos(loc,allele1)]+=difference;
            Epsilon[EpsPos(loc,allele2)]-=difference;
          }
        }
      }
    }
  }
}
