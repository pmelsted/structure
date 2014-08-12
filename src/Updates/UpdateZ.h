extern void UpdateZ (int *Z, float *Q, float *P, int *Geno,float * random);
extern void UpdateZCL (CLDict *clDict,int *Z, float *Q, float *P, int *Geno,
                       float * random);
extern float UpdateZandSingleR (int *Z,  float *Q, float *P, int *Geno,
                                 float *R, float *Mapdistance, int rep, float *Phase,
                                 int *Z1,int *Phasemodel, float *sumIndLikes,
                                 float *indlike_norm);
extern float UpdateZandR (int *Z,  float *Q, float *P, int *Geno,
                           float *R, float *Mapdistance, int rep, float *Phase, int *Z1,
                           int *Phasemodel, float *sumindlike, float *indlike_norm);
