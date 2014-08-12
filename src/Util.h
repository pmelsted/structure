extern void PrintTranslation (int *Translation, int *NumAlleles);
extern void GetNumLocations (struct IND *ind);
extern void Welcome (FILE * file);
extern void Kill ();
extern void CheckParamCombinations();
extern void FreeAll(float *Mapdistance, float *Phase, int *Phasemodel,
                    float *lambda, float *sumlambda,
                    char *Markername, int *Geno, int* PreGeno, int* Recessive,
                    struct IND *Individual,
                    int *Translation, int *NumAlleles, int *Z, int *Z1, float *Q, float *P,
                    float *R, float *sumR, float *varR, float *Epsilon, float *SumEpsilon,
                    float *Fst,
                    float *FstSum, int *NumLociPop, float *PSum, float *QSum,
                    int *AncestDist, float *UsePopProbs, float *LocPrior, float *sumLocPrior,
                    float *Alpha, float *sumAlpha, float *sumIndLikes, float *indLikesNorm,
                    CLDict *clDict);
extern float logsumexp(float a, float b);

