extern void UpdateAlpha (float *Q, float *Alpha, struct IND *Individual,
                         int rep);
extern void UpdateAlphaCL (CLDict *clDict,float *Q, float *Alpha, struct IND *Individual,
                         int rep,int POPFLAGINDS);
extern void UpdateAlphaLocPrior(float *Q, float *Alpha, float *LocPrior,
                                struct IND *Individual);

