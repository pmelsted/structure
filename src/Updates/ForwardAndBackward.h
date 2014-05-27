extern double Forward (int *Z, double *IndividualQ, double *P,
		       int *Geno, double Rec, int ind, double *RTransitProb,
                       double *Mapdistance, double *Phase,int *Phasemodel);


extern void Backward (int *Z, double *IndividualQ,
		      double Rec, int ind, double *Mapdistance,
                      double *RTransitProb, int rep, int *Z1,
		      double *Phase, double *P, int *Geno,int *Phasemodel);

