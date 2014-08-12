#include "../KernelDefs.h"
extern void UpdateQ (int *Geno, int *PreGeno, float *Q, float *P,
                     int *Z, float *Alpha, int rep,
                     struct IND *Individual, float *UsePopProbs,
                     int *Recessive, float *LocPrior,float * randomArr);
extern void UpdateQCL (CLDict *clDict,int *Geno, int *PreGeno, float *Q, float *P,
                     int *Z, float *Alpha, int rep,
                     struct IND *Individual, float *UsePopProbs,
                     int *Recessive, float *LocPrior,float * randomArr);
extern void UpdateQMetroRecombine (int *Geno, float *Q, int *Z, float *P,
                                   float *Alpha, int rep,
                                   struct IND *Individual,
                                   float *Mapdistance, float *R,
                                   float *Phase,int *Phasemodel,
                                   float * randomArr);

