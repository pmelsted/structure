extern float Forward (int *Z, float *IndividualQ, float *P,
                       int *Geno, float Rec, int ind, float *RTransitProb,
                       float *Mapdistance, float *Phase,int *Phasemodel);


extern void Backward (int *Z, float *IndividualQ,
                      float Rec, int ind, float *Mapdistance,
                      float *RTransitProb, int rep, int *Z1,
                      float *Phase, float *P, int *Geno,int *Phasemodel);

