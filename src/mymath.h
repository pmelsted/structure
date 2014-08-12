#include "randGen.h"
extern float Square(float x);
extern float SampleVar(float sumsq,float sum,long num);
extern int PickAnOption(int total,float sum,float cutoffs[]);
extern int PickAnOptionDiscrete(int total,float sum,float Probs[],
                                RndDiscState *randState);
extern float SD(float sumsq, float sum, long num);
extern float LDirichletProb(float prior[],float post[],int length);
extern float LGammaDistProb(float alpha,float beta, float y);
extern float FindAveLogs(float *logmax,float *sum, float lognext,int rep);
extern void RandomOrder(int list[],int length);
extern float Factorial(int n);
extern float mylgamma(float z);
extern float ChiSq(int *list1,int len1,int *list2,int len2,int mincount,
                    int missing,int *df);

