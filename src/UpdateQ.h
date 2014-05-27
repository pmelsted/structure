extern void UpdateQAdmixture (double *Q, int *Z, double *Alpha,
			      struct IND *Individual); 


extern void UpdateQMetro (int *Geno, int *PreGeno, double *Q,
			              double *P, double *Alpha, int rep,
			              struct IND *Individual, int *Recessive);


extern void UpdateQNoAdmix (int *Geno, double *Q, double *P, 
			                struct IND *Individual, double *LocPrior);

extern void UpdateQAdmixture (double *Q, int *Z, double *Alpha,
                              struct IND *Individual);

extern void UpdateQWithPops (int *Geno, double *Q, double *P, int *Z,
                             double *Alpha, int rep, struct IND *Individual,
                             double *UsePopProbs);
                            
extern void UpdateQ (int *Geno, int *PreGeno, double *Q, double *P,
                     int *Z, double *Alpha, int rep,
                     struct IND *Individual, double *UsePopProbs,
                     int *Recessive, double *LocPrior);
                     
extern void UpdateQMetroRecombine (int *Geno, double *Q, int *Z, double *P,
                                   double *Alpha, int rep,
                                   struct IND *Individual,
                                   double *Mapdistance, double *R,
                                   double *Phase,int *Phasemodel);
