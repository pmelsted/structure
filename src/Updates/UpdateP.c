#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"
#include "../Kernels.h"



/*===============================================*/
/*O(MAXPOPS*numalleles + NUMINDS*LINES) */
void
GetNumFromPop (int *NumAFromPop, int *Geno, int *Z, int loc,
               int numalleles,struct IND *Individual)
{
    /*Fill in the number of each allele from each pop */
    int ind, line, pop, allele;
    /* int genpos; */
    int allelevalue;
    int popvalue;

    /* O(MAXPOPS*numalleles) */
    for (pop = 0; pop < MAXPOPS; pop++) {
        for (allele = 0; allele < numalleles; allele++) {
            NumAFromPop[NumAFromPopPos (pop, allele)] = 0;
        }
    }

    /* O(NUMINDS*LINES) */
    if (PFROMPOPFLAGONLY) {     /*this option uses only individuals with POPFLAG=1 to update P*/
        for (ind = 0; ind < NUMINDS; ind++) {
            if (Individual[ind].PopFlag ==
                    1) {    /*individual must have popflag turned on*/
                for (line = 0; line < LINES; line++) {
                    popvalue = Z[ZPos (ind, line, loc)];
                    allelevalue = Geno[GenPos (ind, line, loc)];

                    if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
                        NumAFromPop[NumAFromPopPos (popvalue, allelevalue)]++;
                    }
                }
            }
        }
    } else {       /*standard update--use everybody to update P */
        for (ind = 0; ind < NUMINDS; ind++) {
            for (line = 0; line < LINES; line++) {
                popvalue = Z[ZPos (ind, line, loc)];
                allelevalue = Geno[GenPos (ind, line, loc)];

                if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
                    NumAFromPop[NumAFromPopPos (popvalue, allelevalue)]++;
                }
            }
        }
    }
}

void GetNumFromPopsCL (CLDict *clDict,int *NumAFromPops, int *Geno, int *Z,
                       int *NumAlleles, struct IND *Individual)
{
    cl_int err;
    size_t local;
    size_t *global;
    /* for error handling in kernel */
    int *error;

    int ind;
    int *popflag;


    error = calloc(1,sizeof(int));
    error[0] = 0;
    global = calloc(2,sizeof(size_t));
    global[0] = NUMINDS;
    global[1] = NUMLOCI;

    if(PFROMPOPFLAGONLY) {
        popflag = calloc(NUMINDS,sizeof(int));
        for(ind = 0; ind < NUMINDS; ++ind) {
            popflag[ind] = Individual[ind].PopFlag;
        }

        err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[POPFLAGCL],
                                   CL_TRUE,
                                   0, sizeof(int) * NUMINDS, popflag, 0, NULL, NULL);
        handleCLErr(err, clDict,"UpdateP Error: Failed to write buffer popflags!");
    }

    err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[GENOCL], CL_TRUE,
                               0, sizeof(int) * GENOSIZE, Geno, 0, NULL, NULL);
    handleCLErr(err, clDict,"UpdateP Error: Failed to write buffer geno!");

    err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[ZCL], CL_TRUE,
                               0, sizeof(int) * ZSIZE, Z, 0, NULL, NULL);
    handleCLErr(err, clDict,"UpdateP Error: Failed to write buffer Z!");

    err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[NUMAFROMPOPSCL],
                               CL_TRUE,
                               0, sizeof(int) * NUMLOCI*MAXPOPS*MAXALLELES, NumAFromPops, 0, NULL, NULL);

    handleCLErr(err, clDict,
                "UpdateP Error: Failed to write buffer NumAFromPopsCL!");

    err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[NUMALLELESCL],
                               CL_TRUE,
                               0, sizeof(int) * NUMLOCI, NumAlleles, 0, NULL, NULL);
    handleCLErr(err, clDict,"UpdateP Error: Failed to write buffer NUMALLELES!");
    err = clEnqueueWriteBuffer(clDict->commands, clDict->buffers[ERRORCL], CL_TRUE,
                               0, sizeof(int), error, 0, NULL, NULL);
    handleCLErr(err, clDict,"UpdateP Error: Failed to write error buffer!");



    err = 0;
    err  = clSetKernelArg(clDict->kernels[GetNumFromPopsKernel], 0, sizeof(cl_mem),
                          &(clDict->buffers[ZCL]));
    handleCLErr(err, clDict,"UpdateP Error: Failed to set arg 0!");
    err = clSetKernelArg(clDict->kernels[GetNumFromPopsKernel], 1, sizeof(cl_mem),
                         &(clDict->buffers[GENOCL]));
    handleCLErr(err, clDict,"UpdateP Error: Failed to set arg 1!");
    err = clSetKernelArg(clDict->kernels[GetNumFromPopsKernel], 2, sizeof(cl_mem),
                         &(clDict->buffers[NUMALLELESCL]));
    handleCLErr(err, clDict,"UpdateP Error: Failed to set arg 2!");
    err = clSetKernelArg(clDict->kernels[GetNumFromPopsKernel], 3, sizeof(cl_mem),
                         &(clDict->buffers[POPFLAGCL]));
    handleCLErr(err, clDict,"UpdateP Error: Failed to set arg 3!");

    err = clSetKernelArg(clDict->kernels[GetNumFromPopsKernel], 4, sizeof(cl_mem),
                         &(clDict->buffers[NUMAFROMPOPSCL]));
    handleCLErr(err, clDict,"UpdateP Error: Failed to set arg 4!");

    err = clSetKernelArg(clDict->kernels[GetNumFromPopsKernel], 5, sizeof(cl_mem),
                         &(clDict->buffers[ERRORCL]));
    handleCLErr(err, clDict,"UpdateP Error: Failed to set arg 5!");

    err = clGetKernelWorkGroupInfo(clDict->kernels[GetNumFromPopsKernel],
                                   clDict->device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
    handleCLErr(err, clDict,
                "UpdateP Error: Failed to retrieve kernel work group info!");


    err = clEnqueueNDRangeKernel(clDict->commands,
                                 clDict->kernels[GetNumFromPopsKernel],
                                 2, NULL, global, NULL, 0, NULL, NULL);
    handleCLErr(err, clDict,"UpdateP Error: Failed to execute kernel!");

    err = clFinish(clDict->commands);

    err = clEnqueueReadBuffer(clDict->commands, clDict->buffers[NUMAFROMPOPSCL],
                              CL_TRUE, 0,
                              sizeof(int) * NUMLOCI*MAXPOPS * MAXALLELES, NumAFromPops, 0, NULL, NULL );

    err = clEnqueueReadBuffer(clDict->commands, clDict->buffers[ERRORCL], CL_TRUE,
                              0, sizeof(int), error, 0, NULL, NULL );


    free(global);
    /* some error handling */
    if (error[0] != KERNEL_SUCCESS ) {
        printf("UpdateP Error in Kernel:\n");
        PrintKernelError(error[0]);
        ReleaseCLDict(clDict);
        exit(EXIT_FAILURE);
    }
    free(error);


}

void GetNumFromPops (int *NumAFromPops, int *Geno, int *Z, int *NumAlleles,
                     struct IND *Individual)
{
    int loc;
    int ind, line, pop, allele;
    int allelevalue;
    int popvalue;
    int numalleles;
    int offset;

    for (loc = 0; loc < NUMLOCI; loc++) {
        numalleles = NumAlleles[loc];
        offset = loc*MAXPOPS*MAXALLELES;
        /*Fill in the number of each allele from each pop */
        /* int genpos; */

        /* O(MAXPOPS*numalleles) */
        for (pop = 0; pop < MAXPOPS; pop++) {
            for (allele = 0; allele < numalleles; allele++) {
                NumAFromPops[NumAFromPopPos (pop, allele)+offset] = 0;
            }
        }

        /* O(NUMINDS*LINES) */
        for (ind = 0; ind < NUMINDS; ind++) {
            if (!PFROMPOPFLAGONLY
                    || Individual[ind].PopFlag ==
                    1) {    /*individual must have popflag turned on*/
                for (line = 0; line < LINES; line++) {
                    popvalue = Z[ZPos (ind, line, loc)];
                    allelevalue = Geno[GenPos (ind, line, loc)];

                    if ((allelevalue != MISSING) && (popvalue != UNASSIGNED)) {
                        NumAFromPops[NumAFromPopPos (popvalue, allelevalue) + offset]++;
                    }
                }
            }
        }
    }
}

/*------------------------------------------*/
/*
 * O(NUMLOCI*(MAXPOPS* (max_loc NumAlleles[loc]) + NUMINDS*LINES)) =>
 * O(NUMLOCI*(MAXPOPS*MAXALLELES + NUMINDS*LINES))
 */
void UpdateP (double *P, double *LogP, double *Epsilon, double *Fst,
              int *NumAlleles, int *Geno, int *Z, double *lambda, struct IND *Individual,
              double * randomArr)
/*Simulate new allele frequencies from Dirichlet distribution */
{
    int loc, pop, allele;
    double *Parameters;           /*[MAXALLS] **Parameters of posterior on P */
    /*int *NumAFromPop;             [>[MAXPOPS][MAXALLS] **number of each allele from each pop <]*/

    int *NumAFromPops;/*[NUMLOCI][MAXPOPS][MAXALLS] **number of each allele from each pop at each loc */
    int popsoffset;

    RndDiscState randState[1];

    Parameters = calloc(MAXALLELES, sizeof (double));
    /*NumAFromPop = calloc(MAXPOPS * MAXALLELES, sizeof (int));*/

    NumAFromPops = calloc(NUMLOCI*MAXPOPS * MAXALLELES, sizeof (int));

    if ((Parameters == NULL) || (NumAFromPops == NULL)) {
        printf ("WARNING: unable to allocate array space in UpdateP\n");
        Kill ();
    }

    /*initialize the NumAFromPops array*/
    GetNumFromPops (NumAFromPops,Geno, Z, NumAlleles, Individual);


    /* O(NUMLOCI*(MAXPOPS* (max_loc NumAlleles[loc]) + NUMINDS*LINES)) */

    initRndDiscState(randState,randomArr,MAXALLELES*MAXRANDOM);
    for (loc = 0; loc < NUMLOCI; loc++) {
        popsoffset = loc*MAXPOPS*MAXALLELES;
        /*count number of each allele from each pop */
        /*O(MAXPOPS*NumAlleles[loc] + NUMINDS*LINES) */

        /*
         * Testing:
         */
        /*GetNumFromPop (NumAFromPop, Geno, Z, loc, NumAlleles[loc], Individual);
        for (pop = 0; pop < MAXPOPS; pop++) {
            for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (NumAFromPops[NumAFromPopPos (pop, allele)+popsoffset] !=
                    NumAFromPop[NumAFromPopPos (pop, allele)])
                {
                    printf("NUMAFROMPOPS NOT CORRECT!!!!\n\n\n\n\n\n\n\n\n\n\n\n\n");
                }

            }
        }*/

        /* O(MAXPOPS*NumAlleles[loc])*/
        for (pop = 0; pop < MAXPOPS; pop++) {
            rndDiscStateReset(randState,popsoffset*MAXRANDOM+pop*MAXALLELES*MAXRANDOM);
            for (allele = 0; allele < NumAlleles[loc]; allele++) {
                if (FREQSCORR) {
                    Parameters[allele] = Epsilon[EpsPos (loc, allele)]
                                         *(1.0- Fst[pop])/Fst[pop]
                                         + NumAFromPops[NumAFromPopPos (pop, allele)+popsoffset];
                } else {
                    Parameters[allele] = lambda[pop]
                                         + NumAFromPops[NumAFromPopPos (pop, allele)+popsoffset];
                }
            }
            /*return a value of P simulated from the posterior Di(Parameters) */
            /*O(NumAlleles[loc]) */
            LogRDirichletDisc (Parameters, NumAlleles[loc],
                               P + PPos (loc, pop, 0),
                               LogP +PPos(loc,pop,0),
                               randState);

            /*need to worry about underflow in UpdateEpsilon due to
              allele frequencies being set to zero---hence previously used the
              following hack, however now pass LogP instead

              for (allele=0;allele<NumAlleles[loc];allele++) if
              (P[PPos(loc,pop,allele)]<1E-20) {
              P[PPos(loc,pop,allele)]=1E-20;

              for (pop=0; pop<MAXPOPS; pop++)
              {
              printf(" loc =%d pop= %d fst=%f ",loc,pop,Fst[pop]);
              for (allele=0;allele<NumAlleles[loc];allele++)
              printf (" Epsilon= %.5E P= %.5E Parameters=  %.5E Num= %d",
              Epsilon[EpsPos(loc,allele)],P[PPos(loc,pop,allele)],
              Parameters[allele],NumAFromPop[NumAFromPopPos (pop, allele)]);
              printf("\n");
              }
              } */
        }
    }

    free (Parameters);
    free (NumAFromPops);
}


void UpdatePCL (CLDict *clDict,double *P, double *LogP, double *Epsilon,
                double *Fst,
                int *NumAlleles, int *Geno, int *Z, double *lambda, struct IND *Individual,
                double * randomArr)
/*Simulate new allele frequencies from Dirichlet distribution */
{

    int *NumAFromPops;
    size_t *global;
    /* for error handling in kernel */
    int *error;

    int ind;
    int *popflag;

    NumAFromPops = calloc(NUMLOCI*MAXPOPS * MAXALLELES, sizeof (int));
    error = calloc(2,sizeof(int));
    error[0] = 0;
    error[1] = 0;
    global = calloc(2,sizeof(size_t));
    global[0] = NUMINDS;
    global[1] = NUMLOCI;

    if(PFROMPOPFLAGONLY) {
        popflag = calloc(NUMINDS,sizeof(int));
        for(ind = 0; ind < NUMINDS; ++ind) {
            popflag[ind] = Individual[ind].PopFlag;
        }

        writeBuffer(clDict,popflag,sizeof(int) * NUMINDS,POPFLAGCL,"Popflags");
    }
    /*
     * GetNumFromPops writes
     */

    /* =================================================== */
    writeBuffer(clDict,Geno,sizeof(int) * GENOSIZE,GENOCL,"Geno");
    writeBuffer(clDict,Z,sizeof(int)*ZSIZE,ZCL,"Z");

    writeBuffer(clDict,NumAFromPops,sizeof(int)* NUMLOCI*MAXPOPS*MAXALLELES,
                NUMAFROMPOPSCL,"NumAFromPops");

    writeBuffer(clDict,NumAlleles,sizeof(int) * NUMLOCI,NUMALLELESCL,"NumAlleles");

    writeBuffer(clDict,error,sizeof(int),ERRORCL,"error");

    /* =================================================== */

    /*
     * UpdateP writes
     */

    /* =================================================== */

    readBuffer(clDict,P,sizeof(double) * PSIZE,PCL,"P");
    readBuffer(clDict,LogP,sizeof(double) * PSIZE,LOGPCL,"LogP");

    if (FREQSCORR) {
        writeBuffer(clDict,Fst,sizeof(double) * MAXPOPS,FSTCL,"FST");
        writeBuffer(clDict,Epsilon,sizeof(double) * NUMLOCI*MAXALLELES,EPSILONCL,
                    "EPSILON");
    } else {
        writeBuffer(clDict,lambda,sizeof(double) * MAXPOPS,LAMBDACL,"LAMBDA");
    }

    writeBuffer(clDict,randomArr,
                sizeof(double) * NUMLOCI*MAXALLELES*MAXPOPS*MAXRANDOM,RANDCL,"randomArr");


    /* =================================================== */



    runKernel(clDict,GetNumFromPopsKernel,2,global,"GetNumFromPops");


    global[0] = NUMLOCI;
    global[1] = MAXPOPS;

    runKernel(clDict,UpdatePKernel,2,global,"UpdateP");

    readBuffer(clDict,P,sizeof(double) * PSIZE,PCL,"P");
    readBuffer(clDict,LogP,sizeof(double) * PSIZE,LOGPCL,"LogP");
    readBuffer(clDict,error,sizeof(int),ERRORCL,"Error");



    free(global);
    /* some error handling */
    if (error[0] != KERNEL_SUCCESS ) {
        printf("UpdateP Error in Kernel:\n");
        PrintKernelError(error[0]);
        printf("%d\n",error[1]);
        ReleaseCLDict(clDict);
        exit(EXIT_FAILURE);
    }
    free(error);
    /*GetNumFromPops (NumAFromPops,Geno, Z, NumAlleles, Individual);
    for(loc = 0; loc < NUMLOCI*MAXALLELES*MAXPOPS; ++loc){
        if(NumAFromPopsCL[loc] != NumAFromPops[loc]){
            printf("NumAFromPopsCL not correct!\n");
            printf("loc %d, size %d\n", loc, NUMLOCI*MAXALLELES*MAXRANDOM);
            printf("%d, %d\n",NumAFromPopsCL[loc],NumAFromPops[loc]);
            ReleaseCLDict(clDict);
            exit(EXIT_FAILURE);
        }
    }*/


    free (NumAFromPops);
}
