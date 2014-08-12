#include "KernelDefs.h"
extern float CalcLikeInd(int *Geno, int *PreGeno, float *Q, float *P,
                   int ind, int *Recessive);

extern float CalcLikeIndDiff(int *Geno, int *PreGeno, float *QPlus, float *QMinus,
                       float *P,
                       int ind, int *Recessive);

extern float CalcLikeIndDiffCL(int *Geno,  float *TestQ, float *Q,
                       float *P, int ind);
extern float CalcLikeIndCL(int *Geno, float *Q, float *P, int ind);
extern float CalcLike (int *Geno, int *PreGeno, float *Q, float *P, int *Recessive,
                 float *sumindlike, float *indlike_norm);
extern void CalcLogdiffsCL(CLDict *clDict,int *Geno,float *TestQ, float *Q, float *P, float *logdiffs);

