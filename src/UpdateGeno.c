#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ran.h"
#include "mymath.h"
#include "structure.h"
#include "output.h"

/*------------------------------------------*/
void UpdateGeno (int *PreGeno, int *Geno, double *P, int *Z,
                 int *Recessive, int *NumAlleles, double *Q)
    /* 
     * this function updates the imputed genotypes when the genotypes are
     * ambiguous due to recessive markers or inbreeding.
     */
{

  int ind, loc, allele, line, dom, toggle,rejectioncount,allelecount,notmissingcount;
  int *AlleleUsed=NULL, *AllelePresent=NULL;
  double *AlleleProbs, Sum;
  static int RejectionThreshold = 1000000;
  /*  int pop;
   *  double temp, Sum1, Sum2; */

  if (LINES == 2) {
    AlleleProbs = calloc (4, sizeof (double));
  } else {
    AlleleProbs = calloc (MAXALLELES, sizeof (double));
    AlleleUsed = calloc (MAXALLELES, sizeof (int));
    AllelePresent = calloc (MAXALLELES, sizeof (int));
  }

  if (LINES == 2) {  /* special update for diploids; general case by rejection sampling */
    for (ind = 0; ind < NUMINDS; ind++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        if (PreGeno[GenPos (ind, 0, loc)] != MISSING
            && PreGeno[GenPos (ind, 1, loc)] != MISSING) {
          if (PreGeno[GenPos (ind, 0, loc)] ==
              PreGeno[GenPos (ind, 1, loc)]) {
            for (dom = 0; dom < 4; dom++) {
              AlleleProbs[dom] = 0.0;
            }

            if (RECESSIVEALLELES
                && Recessive[loc] != MISSING    /* bug fixed 05072007 */
                && PreGeno[GenPos (ind, 0, loc)] != Recessive[loc]) {
              AlleleProbs[0] =
                  P[PPos(loc, Z[ZPos(ind, 0, loc)], Recessive[loc])] *
                  P[PPos(loc, Z[ZPos(ind, 1, loc)], PreGeno[GenPos(ind, 0, loc)])];
              AlleleProbs[1] =
                  P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                  P[PPos(loc, Z[ZPos (ind, 1, loc)], Recessive[loc])];
            }

            AlleleProbs[2] =
                P[PPos(loc, Z[ZPos (ind, 0, loc)], PreGeno[GenPos (ind, 0, loc)])] *
                P[PPos(loc, Z[ZPos (ind, 1, loc)], PreGeno[GenPos (ind, 1, loc)])];

            Sum = AlleleProbs[0] + AlleleProbs[1] + AlleleProbs[2];
            dom = PickAnOption (3, Sum, AlleleProbs);

            if (dom == 0) {
              Geno[GenPos (ind, 0, loc)] = Recessive[loc];
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else if (dom == 1) {
              Geno[GenPos (ind, 1, loc)] = Recessive[loc];
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else if (dom == 2) {
              Geno[GenPos (ind, 1, loc)] = PreGeno[GenPos (ind, 0, loc)];
              Geno[GenPos (ind, 0, loc)] = PreGeno[GenPos (ind, 0, loc)];
            } else {
              printf ("surely some mistake in UpdateGeno\n");
              Kill(); /* modify by William - kill is not a standard ANSI C function */
            }
          }
        }
      }
    }
  } else { /* general n-ploid case.  Select genotype by rejection sampling */
    for (loc = 0; loc < NUMLOCI; loc++) {
      for (ind = 0; ind < NUMINDS; ind++) {
        if (Recessive[loc]==NOTAMBIGUOUS) {
          for (line=0;line<LINES;line++) {
            Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
          }
        } else {
          allelecount=0;
          notmissingcount=0;
          
          for (allele = 0; allele < NumAlleles[loc]; allele++) {
            AllelePresent[allele] = 0;
          }
          
          for (line = 0; line < LINES; line++) {
            if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
              notmissingcount+=1;
              if (AllelePresent[PreGeno[GenPos (ind, line, loc)]]==0) {
                AllelePresent[PreGeno[GenPos (ind, line, loc)]] = 1;
                allelecount+=1;
              }
            }
          }
          
          if (allelecount==notmissingcount) {  /* if number of alleles equal to number of slots then nothing to do */
            for (line=0;line<LINES;line++) {
              Geno[GenPos(ind,line,loc)]=PreGeno[GenPos(ind,line,loc)];
            }
          } else {
            toggle = 1;
            rejectioncount=0;
            while (toggle && rejectioncount < RejectionThreshold) {
              rejectioncount+=1;
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                AlleleUsed[allele] = 0;
              }
              
              for (line = 0; line < LINES; line++) {
                if (PreGeno[GenPos (ind, line, loc)] != MISSING) {
                  for (allele = 0; allele < NumAlleles[loc]; allele++) {
                    if (AllelePresent[allele] || allele == Recessive[loc]) {
                      AlleleProbs[allele] =
                          P[PPos (loc, Z[ZPos (ind, line, loc)], allele)];
                    } else {
                      AlleleProbs[allele] = 0.0;
                    }
                  }
                  Sum = 0.0;
                  for (allele = 0; allele < NumAlleles[loc]; allele++) {
                    Sum += AlleleProbs[allele];
                  }
                  dom = PickAnOption (NumAlleles[loc], Sum, AlleleProbs);
                  Geno[GenPos (ind, line, loc)] = dom;
                  AlleleUsed[dom] = 1;
                }
              }
              
              toggle = 0;
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (AlleleUsed[allele] == 0 && AllelePresent[allele] == 1) {
                  toggle = 1;
                }
              }
            }
            
            if (toggle==1) {
              /* rejection method failed, set a lower rejection threshold for future iterations */
              if (RejectionThreshold > 100) {
                RejectionThreshold /= 2;
                if (RejectionThreshold <100) {
                  RejectionThreshold = 100;
                }
              } 
              /*
                printf("sorry, STRUCTURE has tried to simulate alleles for individual %d at locus %d 1,000,000 times by a rejection method and failed", ind+1, loc+1);
                Kill();
              */
              line=0;

              /*  first fill in one copy of each present allele. */
              for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (AllelePresent[allele] == 1) {
                  Geno[GenPos(ind,line,loc)]= allele;
                  line++;
                }
              }

              /* now fill in remaining nonmissing slots. */
              for (line=allelecount; line<notmissingcount; line++) {
                for (allele = 0; allele < NumAlleles[loc]; allele++) {
                  if (AllelePresent[allele] || allele == Recessive[loc]) {
                    AlleleProbs[allele] = P[PPos (loc, Z[ZPos (ind, line, loc)], allele)];
                  } else {
                    AlleleProbs[allele] = 0.0;
                  }
                }

                Sum = 0.0;
                for (allele = 0; allele < NumAlleles[loc]; allele++) {
                  Sum += AlleleProbs[allele];
                }
                dom = PickAnOption (NumAlleles[loc], Sum, AlleleProbs);
                Geno[GenPos (ind, line, loc)] = dom;
                AlleleUsed[dom] = 1;
              }

              /* now fill in missing slots */
              for (line=notmissingcount;line<LINES;line++) {
                Geno[GenPos (ind, line, loc)]=MISSING;
              }
            }
          }
        }
      }
    }
  }
  free (AlleleProbs);
  if (LINES != 2) {
    free (AlleleUsed);
    free (AllelePresent);
  }
}
