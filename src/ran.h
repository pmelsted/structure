#include "randGen.h"
extern void Randomize(int RANDOMIZE, int *seed);
extern float RandomReal(float low, float high);
extern int RandomInteger(int low, int high);
extern float rnd();
extern float RGamma(float n,float lambda);
extern void RDirichlet(const float * a, const int k, float * b);
extern long RPoisson(float mu);
extern float RExpon(float av);
extern float RNormal(float mu,float sd) ;
extern float fsign( float num, float sign );
extern float sexpo(void);
extern float snorm();
extern float genexp(float av);
extern long ignpoi(float mean);
extern long ignuin(int low, int high);
extern float genunf(float low, float high);
extern long   Binomial(int n, float p);
extern long   Binomial1(int n, float p);
extern float BinoProb(int n, float p,int i);
extern void LogRDirichlet (const float *a, const int k, float *b,float *c);
extern float numToRange(float low, float high, float num);
extern float rndDisc(RndDiscState *randState);
extern float RGammaDisc(float n,float lambda,RndDiscState *randState);
extern void RDirichletDisc(const float * a, const int k, float * b, int offset,
                           RndDiscState *randState);
extern float LogRGammaDisc(float n, float lambda, RndDiscState *randState);
extern void LogRDirichletDisc (const float *a, const int k, float *b,
                                RndDiscState *randState);

