extern void PrintTranslation (int *Translation, int *NumAlleles);
extern void GetNumLocations (struct IND *ind);
extern void Welcome (FILE * file);
extern void Kill ();
extern void CheckParamCombinations();
extern void FreeAll(double *Mapdistance, double *Phase, int *Phasemodel, double *lambda, double *sumlambda,
	     char *Markername, int *Geno, int* PreGeno, int* Recessive, struct IND *Individual,
	     int *Translation, int *NumAlleles, int *Z, int *Z1, double *Q, double *P, double *LogP,
	     double *R, double *sumR, double *varR, double *Epsilon, double *SumEpsilon, double *Fst,
	     double *FstSum, int *NumLociPop, double *PSum, double *QSum,
	     int *AncestDist, double *UsePopProbs, double *LocPrior, double *sumLocPrior,
	     double *Alpha, double *sumAlpha, double *sumIndLikes, double *indLikesNorm);
extern double logsumexp(double a, double b);
