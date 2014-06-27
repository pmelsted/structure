#include "KernelDefs.h"
extern double CalcLikeInd(int *Geno, int *PreGeno, double *Q, double *P,
                   int ind, int *Recessive);

extern double CalcLikeIndDiff(int *Geno, int *PreGeno, double *QPlus, double *QMinus,
                       double *P,
                       int ind, int *Recessive);

extern double CalcLikeIndDiffCL(int *Geno,  double *TestQ, double *Q,
                       double *P, int ind);
extern double CalcLikeIndCL(int *Geno, double *Q, double *P, int ind);
extern double CalcLike (int *Geno, int *PreGeno, double *Q, double *P, int *Recessive,
                 double *sumindlike, double *indlike_norm);
extern void CalcLogdiffsCL(CLDict *clDict,int *Geno,double *TestQ, double *Q, double *P, double *logdiffs);

