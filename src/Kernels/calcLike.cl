/*
    This version adds multiple elements per thread sequentially.  This reduces the overall
    cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
    (Brent's Theorem optimization)
*/
__kernel void reduce(__global double *g_idata, __global double *g_odata, unsigned int n, __local volatile double *sdata)
{
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);
    unsigned int gridSize = blockSize*2*get_num_groups(0);
    sdata[tid] = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridDim).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        sdata[tid] += g_idata[i];
        // ensure we don't read out of bounds -- this is optimized away for powerOf2 sized arrays
        if (i + blockSize < n){
            sdata[tid] += g_idata[i+blockSize];
        }

        i += gridSize;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // do reduction in shared mem
    if (blockSize >= 512) {
        if (tid < 256) { sdata[tid] += sdata[tid + 256]; }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } barrier(CLK_LOCAL_MEM_FENCE); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } barrier(CLK_LOCAL_MEM_FENCE); }

    if (tid < 32)
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; }
    }

    // write result for this block to global mem
    if (tid == 0) g_odata[get_group_id(0)] = sdata[0];
}


__kernel void reduceLogDiffs(__global double *logterms,
                             __global double *logdiffs,
                             __local  double *scratch,
                             const int ind)
{
    int loc = get_global_id(1);
    double logdiff = 0.0;
    while( loc < NUMLOCI){
        double elem = logterms[ind*NUMLOCI+loc];
        logdiff += elem;
        loc += get_global_size(1);
    }

    int localLoc = get_local_id(1);
    scratch[localLoc] = logdiff;
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int offset = get_local_size(1) /2; offset > 0; offset >>= 1){
        if(localLoc < offset){
            scratch[localLoc] += scratch[localLoc + offset];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(localLoc == 0){
        logdiffs[get_group_id(1)] = scratch[0];
    }
}


double mapLogDiffsFunc(__global double *Q, __global double *TestQ,
                       __global double *P, __global int *Geno,
                       int ind, int loc)
{
    int allele, line, pop;
    double termP;
    double termM;
    double logterm = 0.0;

    for (line = 0; line < LINES; line++) {
        allele = Geno[GenPos (ind, line, loc)];
        if (allele != MISSING) {
            termP = 0.0;
            termM = 0.0;
            for (pop = 0; pop < MAXPOPS; pop++) {
                termP += TestQ[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                termM += Q[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
            }
            //logterms[ind*NUMLOCI + loc] += log(termP) - log(termM);
            logterm += log(termP) - log(termM);
        }
    }
    return logterm;
}



/* ===== UNFINISHED ======= */
__kernel void mapReduceLogDiffs(__global double *Q,
                                __global double *TestQ,
                                __global double *P,
                                __global int *Geno,
                                __global double *logdiffs,
                                __global double *results,
                                __local  double *scratch)
{
    int ind = get_global_id(0);
    int loc = get_global_id(1);
    int globalsize = get_global_size(1);
    double logdiff = 0.0;
    /* ======= Map to functions, and store in local memory ======== */
    while( loc < NUMLOCI){
        double elem = mapLogDiffsFunc(Q,TestQ,P,Geno,ind,loc)
        logdiff += elem;
        loc += globalsize;
    }

    int localLoc = get_local_id(1);
    scratch[localLoc] = logdiff;
    barrier(CLK_LOCAL_MEM_FENCE);
    /* ======== Reduce local memory ==========*/
    for(int offset = get_local_size(1) /2; offset > 0; offset >>= 1){
        if(localLoc < offset){
            scratch[localLoc] += scratch[localLoc + offset];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    int group = get_group_id(1);
    if(localLoc == 0){
        results[group] = scratch[0];
        barrier(CLK_GLOBAL_MEM_FENCE);
        if (group == 0){
            for(group =0; group < globalsize; group++){
                logdiffs[ind] += results[group];
            }
        }
    }
}


__kernel void mapLogDiffs(__global double *logterms, __global double *Q,
                          __global double *TestQ, __global double *P,
                          __global int *Geno, __global int *error)
{

    int ind = get_global_id(0);
    int loc = get_global_id(1);

    logterms[ind*NUMLOCI + loc] = mapLogDiffsFunc(Q,TestQ,P,Geno,ind,loc);

}




