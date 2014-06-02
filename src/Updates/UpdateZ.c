#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ran.h"
#include "../mymath.h"
#include "../structure.h"
#include "../Kernels.h"

#include "ForwardAndBackward.h"

/*-----------------------------------------*/
/*O*(NUMINDS*LINES*NUMLOCI*MAXPOPS)*/
void UpdateZ (int *Z,  double *Q, double *P, int *Geno,double * randomArr)
    /*update Z: population origin of each allele */
{
  int ind, line, loc, pop;
  double *Cutoffs /*[MAXPOPS] */ ;
  /*Cutoffs contains unnormalized probabilities of
    an allele coming from each population */
  double sum=0.0;
  int allele;

  Cutoffs = calloc (MAXPOPS, sizeof (double));
  if (Cutoffs == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZ\n");
    Kill ();
  }

  /*O*(NUMINDS*LINES*NUMLOCI*MAXPOPS)*/
  for (ind = 0; ind < NUMINDS; ind++) {  /*go through all alleles in sample */
    for (line = 0; line < LINES; line++) {
      for (loc = 0; loc < NUMLOCI; loc++) {
        
        allele = Geno[GenPos (ind, line, loc)];

        if (allele == MISSING) {
          Z[ZPos (ind, line, loc)] = UNASSIGNED;
        } else {
          sum = 0.0f;
          for (pop = 0; pop < MAXPOPS; pop++) {
            Cutoffs[pop] = Q[QPos (ind, pop)] * P[PPos (loc, pop, allele)];
            sum += Cutoffs[pop];
          }

          Z[ZPos (ind, line, loc)] = PickAnOptionDiscrete (MAXPOPS, sum, Cutoffs,randomArr[ind*NUMLOCI+loc]);
        }
        
      }
    }
  }

  free (Cutoffs);
}

/*-----------------------------------------*/
/*O*(NUMINDS*LINES*NUMLOCI*MAXPOPS)*/
void UpdateZCL (CLDict *clDict,int *Z,  double *Q, double *P, int *Geno,double * randomArr)
    /*update Z: population origin of each allele */
{
    cl_mem qCL,pCL,genoCL,randCL,zCL;
    cl_int err;
    int i;
    size_t local;                      
    size_t *global;



    global = calloc(2,sizeof(size_t));
    global[0] = NUMINDS;
    global[1] = NUMLOCI;


    qCL = clCreateBuffer(clDict->context,  CL_MEM_READ_ONLY,  sizeof(double)*QSIZE,NULL, &err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed create buffer Q!\n");
        printCLErr(err);
        exit(1);
    }


    pCL = clCreateBuffer(clDict->context,  CL_MEM_READ_ONLY,  sizeof(double)*PSIZE,NULL, &err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed create buffer P!\n");
        printCLErr(err);
        exit(1);
    }
    
    genoCL = clCreateBuffer(clDict->context,  CL_MEM_READ_ONLY,  sizeof(int)*GENOSIZE,NULL, &err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed create buffer Geno!\n");
        printCLErr(err);
        exit(1);
    }

    randCL = clCreateBuffer(clDict->context,  CL_MEM_READ_ONLY,  sizeof(double)*RANDSIZE,NULL, &err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed create buffer Rand!\n");
        printCLErr(err);
        exit(1);
    }

    zCL = clCreateBuffer(clDict->context,  CL_MEM_WRITE_ONLY,  sizeof(int)*ZSIZE,NULL, &err);
    if (err != CL_SUCCESS) {
        printf("Error: Failed create buffer Z!\n");
        printCLErr(err);
        exit(1);
    }

    err = 0;
    err = clEnqueueWriteBuffer(clDict->commands, qCL, CL_TRUE, 0, sizeof(double) * QSIZE, Q, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to write buffer Q!\n");
        printCLErr(err);
        exit(1);
    }
    err = clEnqueueWriteBuffer(clDict->commands, pCL, CL_TRUE, 0, sizeof(double) * PSIZE, P, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to write buffer P!\n");
        printCLErr(err);
        exit(1);
    }
    err = clEnqueueWriteBuffer(clDict->commands, genoCL, CL_TRUE, 0, sizeof(int) * GENOSIZE, Geno, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to write buffer geno!\n");
        printCLErr(err);
        exit(1);
    }
    err = clEnqueueWriteBuffer(clDict->commands, randCL, CL_TRUE, 0, sizeof(double) * RANDSIZE, randomArr, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to write buffer rand!\n");
        printCLErr(err);
        exit(1);
    }


    err = 0;
    err  = clSetKernelArg(clDict->kernels[UpdateZKernel], 0, sizeof(cl_mem), &qCL);
    err |= clSetKernelArg(clDict->kernels[UpdateZKernel], 1, sizeof(cl_mem), &pCL);
    err |= clSetKernelArg(clDict->kernels[UpdateZKernel], 2, sizeof(cl_mem), &genoCL);
    err |= clSetKernelArg(clDict->kernels[UpdateZKernel], 3, sizeof(cl_mem), &randCL);
    err |= clSetKernelArg(clDict->kernels[UpdateZKernel], 4, sizeof(cl_mem), &zCL);

    if (err != CL_SUCCESS) {
        printf("Error: Failed set args! %d\n", err);
        exit(1);
    }

    err = clGetKernelWorkGroupInfo(clDict->kernels[UpdateZKernel], clDict->device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
    if (err != CL_SUCCESS)
    {
        printf("Error: Failed to retrieve kernel work group info! %d\n", err);
        exit(1);
    }
    /*err = clEnqueueNDRangeKernel(clDict->commands, kernel, 2, NULL, &global, &local, 0, NULL, NULL);*/
    err = clEnqueueNDRangeKernel(clDict->commands, clDict->kernels[UpdateZKernel], 2, NULL, global, NULL, 0, NULL, NULL);
    if (err)
    {
        printf("Error: Failed to execute kernel!\n");
        printCLErr(err);
        exit(EXIT_FAILURE);
    }

    err = clFinish(clDict->commands);

    err = clEnqueueReadBuffer(clDict->commands, zCL, CL_TRUE, 0, sizeof(int) * ZSIZE, Z, 0, NULL, NULL );  
    free(global);
    clReleaseMemObject(qCL);
    clReleaseMemObject(pCL);
    clReleaseMemObject(zCL);
    clReleaseMemObject(genoCL);
    clReleaseMemObject(randCL);
}






/*----------------------------------------*/
double
UpdateZandSingleR (int *Z,  double *Q, double *P, int *Geno,
                   double *R, double *Mapdistance, int rep, double *Phase,
                   int *Z1,int *Phasemodel, double *sumIndLikes,
                   double *indlike_norm)
    /* updates Z and R, assuming that the data is phased */
{
  long ind, pop;
  double *RTransitProb;
  double *IndividualQ;
  double trialR, logtrialR,currentloglikelihood, trialloglikelihood, indlike;
  /* long loc; */
  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  IndividualQ = calloc (MAXPOPS, sizeof (double));
  /*this form ensures compatibility with UpdateQMetroRecombine */

  if (RTransitProb == NULL || IndividualQ == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZandSingleR\n");
    Kill ();
  }
  
  currentloglikelihood = 0.0;
  trialloglikelihood = 0.0;
  logtrialR = RNormal(log(R[0])/2.30259,LOG10RPROPSD);
  if (logtrialR<LOG10RMIN) {
    logtrialR=2*LOG10RMIN-logtrialR;
  }

  if (logtrialR>LOG10RMAX) {
    logtrialR=2*LOG10RMAX-logtrialR;
  }
  
  trialR=exp(2.30259*logtrialR);
  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      IndividualQ[pop] = Q[QPos (ind, pop)];
    }
    indlike = Forward(Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                      Mapdistance, Phase,Phasemodel);
    currentloglikelihood += indlike;
    if (sumIndLikes!=NULL) {
      sumIndLikes[ind] += exp(indlike-indlike_norm[ind]);
    }

    Backward(Z,  IndividualQ, R[ind], ind, Mapdistance,
             RTransitProb, rep, Z1, Phase, P, Geno,Phasemodel);

    trialloglikelihood += Forward (Z, IndividualQ, P, Geno, trialR, ind, RTransitProb,
                                   Mapdistance, Phase,Phasemodel);

  }
  /*printf("% .3E % 3E % .8E % .8E % .5E\n",trialR,R[0],trialloglikelihood,currentloglikelihood,exp (trialloglikelihood - currentloglikelihood)); */
  if (RandomReal (0.0, 1.0) < exp (trialloglikelihood - currentloglikelihood)) {  /*Accept */
    R[0] = trialR;
    /*currentloglikelihood=trialloglikelihood;  commented out by JKP--see email from Daniel 9 Dec 02*/
  }
  
  for (ind = 0; ind < NUMINDS; ind++) {
    R[ind] = R[0];
  }
  
  free (RTransitProb);
  free (IndividualQ);

  return currentloglikelihood;
}



/*----------------------------------------*/
double
UpdateZandR (int *Z,  double *Q, double *P, int *Geno,
             double *R, double *Mapdistance, int rep, double *Phase, int *Z1,int *Phasemodel, double *sumindlike, double *indlike_norm)
    /* updates Z and R, assuming that the data is phased */
{
  long ind, pop;
  double *RTransitProb;
  double *IndividualQ;
  double trialR,logtrialR, currentloglikelihood, trialloglikelihood,sumlikelihood;
  /*  long loc; */


  if (PHASED) {
    RTransitProb = calloc (MAXPOPS * NUMLOCI * LINES, sizeof (double));
  } else {
    RTransitProb = calloc (MAXPOPS * MAXPOPS * NUMLOCI, sizeof (double));
  }
  
  IndividualQ = calloc (MAXPOPS, sizeof (double));
  /*this form ensures compatibility with UpdateQMetroRecombine */

  if (RTransitProb == NULL || IndividualQ == NULL) {
    printf ("WARNING: unable to allocate array space in UpdateZandR\n");
    Kill ();
  }
  
  sumlikelihood=0.0;
  for (ind = 0; ind < NUMINDS; ind++) {
    for (pop = 0; pop < MAXPOPS; pop++) {
      IndividualQ[pop] = Q[QPos (ind, pop)];
    }

    currentloglikelihood = Forward (Z, IndividualQ, P, Geno, R[ind], ind, RTransitProb,
                                    Mapdistance, Phase,Phasemodel);
    Backward (Z,  IndividualQ, R[ind], ind, Mapdistance, RTransitProb,
              rep, Z1, Phase, P, Geno,Phasemodel);
    
    logtrialR = RNormal(log(R[ind])/2.30259,LOG10RPROPSD);
    if (logtrialR<LOG10RMIN) {
      logtrialR=2*LOG10RMIN-logtrialR;
    }
    if (logtrialR>LOG10RMAX) {
      logtrialR=2*LOG10RMAX-logtrialR;
    }
    trialR=exp(2.30259*logtrialR);
    trialloglikelihood = Forward (Z, IndividualQ, P, Geno, trialR, ind, RTransitProb, Mapdistance, Phase,Phasemodel);
    /*printf("% .3E % 3E % .8E % .8E % .5E\n",trialR,R[ind],trialloglikelihood,currentloglikelihood,exp (trialloglikelihood - currentloglikelihood)); */
    if (RandomReal (0.0, 1.0) < exp (trialloglikelihood - currentloglikelihood)) {        /*Accept */
      R[ind] = trialR;
      sumlikelihood+=trialloglikelihood;
      if (sumindlike!=NULL) {
        sumindlike[ind] += exp(trialloglikelihood-indlike_norm[ind]);
      }
    } else {
      sumlikelihood+=currentloglikelihood;
      if (sumindlike!=NULL) {
        sumindlike[ind] += exp(currentloglikelihood-indlike_norm[ind]);
      }
    }
  }
  free (RTransitProb);
  free (IndividualQ);

  return sumlikelihood;
}
