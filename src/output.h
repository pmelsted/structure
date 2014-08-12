/* #include "structure.h" */
void
DataCollection (int *Geno,int *PreGeno, float *Q, float *QSum, int *Z,
                int *Z1,  float *P, float *PSum,
                float *Fst, float *FstSum, int *NumAlleles,
                int *AncestDist, float *Alpha, float *sumAlpha,
                float *sumR, float *varR, float *like,
                float *sumlikes, float *sumsqlikes, float *R,
                float *Epsilon, float *SumEpsilon, float recomblikelihood,
                float *lambda, float *sumlambda, int *Recessive,
                float *PopPrior, float *PopPriorSum, int PopPriorLen,
                float *sumindlikes, float *indlikes_norm, int rep);
void
DataCollectionCL (CLDict *clDict,int *Geno,int *PreGeno, float *Q, float *QSum, int *Z,
                int *Z1,  float *P, float *PSum,
                float *Fst, float *FstSum, int *NumAlleles,
                int *AncestDist, float *Alpha, float *sumAlpha,
                float *sumR, float *varR, float *like,
                float *sumlikes, float *sumsqlikes, float *R,
                float *Epsilon, float *SumEpsilon, float recomblikelihood,
                float *lambda, float *sumlambda, int *Recessive,
                float *PopPrior, float *PopPriorSum, int PopPriorLen,
                float *sumindlikes, float *indlikes_norm, int rep);

void PrintLike (float like, int rep, int *Geno, int *PreGeno, float *Q,
                float *P,
                float recomblikelihood);
void PrintUpdate (int rep, int *Geno, int *PreGeno, float *Alpha,
                  float *Correls,
                  float *P, float *Q, float like, float sumlikes,
                  float sumsqlikes, int *NumAlleles, float *R, float *lambda,
                  struct IND *Individual,  float recomblikelihood, int *Recessive,
                  float *PopPrior, int PopPriorLen);
void
OutPutResults (int *Geno, int rep, int savefreq, struct IND *Individual,
               float *PSum, float *QSum,  float *FstSum,
               int *AncestDist, float *UsePopProbs,
               float sumlikes, float sumsqlikes, float *sumAlpha,
               float *sumR, float *varR,
               int *NumAlleles, int *Translation, int final,
               char *Markername, float *R, float *SumEpsilon,
               float *lambda, float *sumlambda, float *PopPriorSum,
               int PopPriorLen, float *sumindlikes, float *indlikes_norm,
               int argc, char *argv[]);

