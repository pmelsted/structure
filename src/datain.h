#include "structure.h"

void ReadInputFile (int *Geno, float *Mapdistances, char *Markernames,
                    struct IND *Individual,float *Phase,int *Recessive);
void CountAlleles (int *Geno, int *NumAlleles, int *Translation,
                   int *Recessive);
