#include "randGen.h"
extern void Randomize(int RANDOMIZE, int *seed);
extern double RandomReal(double low, double high);
extern int RandomInteger(int low, int high);
extern double rnd();
extern double RGamma(double n,double lambda);
extern void RDirichlet(const double * a, const int k, double * b);
extern long RPoisson(double mu);
extern double RExpon(double av);
extern double RNormal(double mu,double sd) ;
extern double fsign( double num, double sign );
extern double sexpo(void);
extern double snorm();
extern double genexp(double av);
extern long ignpoi(double mean);
extern long ignuin(int low, int high);
extern double genunf(double low, double high);
extern long   Binomial(int n, double p);
extern long   Binomial1(int n, double p);
extern double BinoProb(int n, double p,int i);
extern void LogRDirichlet (const double *a, const int k, double *b,double *c);
extern double numToRange(double low, double high, double num);
extern double rndDisc(RndDiscState *randState);
extern double RGammaDisc(double n,double lambda,RndDiscState *randState);
extern void RDirichletDisc(const double * a, const int k, double * b, int offset,
                           RndDiscState *randState);
extern double LogRGammaDisc(double n, double lambda, RndDiscState *randState);
extern void LogRDirichletDisc (const double *a, const int k, double *b,
                                RndDiscState *randState);

