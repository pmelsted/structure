#include "../KernelDefs.h"
extern void UpdateQ (CLDict *clDict,int *Geno, int *PreGeno, double *Q, double *P,
                     int *Z, double *Alpha, int rep,
                     struct IND *Individual, double *UsePopProbs,
                     int *Recessive, double *LocPrior,double * randomArr);
extern void UpdateQMetroRecombine (int *Geno, double *Q, int *Z, double *P,
                                   double *Alpha, int rep,
                                   struct IND *Individual,
                                   double *Mapdistance, double *R,
                                   double *Phase,int *Phasemodel,
                                   double * randomArr);
