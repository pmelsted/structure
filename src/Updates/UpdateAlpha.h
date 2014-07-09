extern void UpdateAlpha (double *Q, double *Alpha, struct IND *Individual,
                         int rep);
extern void UpdateAlphaCL (CLDict *clDict,double *Q, double *Alpha, struct IND *Individual,
                         int rep,int POPFLAGINDS);
extern void UpdateAlphaLocPrior(double *Q, double *Alpha, double *LocPrior,
                                struct IND *Individual);

