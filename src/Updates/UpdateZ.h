extern void UpdateZ (int *Z, double *Q, double *P, int *Geno,double * random);
extern void UpdateZCL (CLDict *clDict,int *Z, double *Q, double *P, int *Geno,double * random);
extern double UpdateZandSingleR (int *Z,  double *Q, double *P, int *Geno,
                                 double *R, double *Mapdistance, int rep, double *Phase,
                                 int *Z1,int *Phasemodel, double *sumIndLikes,
                                 double *indlike_norm);
extern double UpdateZandR (int *Z,  double *Q, double *P, int *Geno,
                           double *R, double *Mapdistance, int rep, double *Phase, int *Z1,
			   int *Phasemodel, double *sumindlike, double *indlike_norm);
