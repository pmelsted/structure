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



double mapLogDiffsFunc(__global double *Q, __global double *TestQ,
                       __global double *P, __global int *Geno,
                       int ind, int loc)
{
    int allele, line, pop;
    double termP;
    double termM;
    if (ind < NUMINDS && loc < NUMLOCI){
        double runningtotalP = 1.0, runningtotalM = 1.0;
        double logtermP = 0.0, logtermM = 0.0;
        for (line = 0; line < LINES; line++) {
            allele = Geno[GenPos (ind, line, loc)];
            if (allele != MISSING) {
                termP = 0.0;
                termM = 0.0;
                for (pop = 0; pop < MAXPOPS; pop++) {
                    termP += TestQ[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                    termM += Q[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                }

                //TODO: Evaluate underflow safe vs nonsafe
                // safe version, should not underflow
                if (termP > SQUNDERFLO) {
                    runningtotalP *= termP;
                } else {
                    runningtotalP *= SQUNDERFLO;
                }
                if (runningtotalP < SQUNDERFLO){
                    logtermP += log(runningtotalP);
                    runningtotalP = 1.0;
                }

                if (termM > SQUNDERFLO) {
                    runningtotalM *= termM;
                } else {
                    runningtotalM *= SQUNDERFLO;
                }
                if (runningtotalM < SQUNDERFLO){
                    logtermM += log(runningtotalM);
                    runningtotalM = 1.0;
                }
                //Might underflow?
                /*logtermP += log(termP);
                logtermM += log(termM);*/
            }
        }
        logtermP += log(runningtotalP);
        logtermM += log(runningtotalM);

        return logtermP - logtermM;
    }
    return 0.0;
}



__kernel void mapReduceLogDiffs(__global double *Q,
                                __global double *TestQ,
                                __global double *P,
                                __global int *Geno,
                                __global double *logdiffs,
                                __global double *results,
                                __local  double *scratch)
{
    int loc = get_global_id(0);
    int ind = get_global_id(1);
    int numgroups = get_num_groups(0);
    /* idempotent */
    double logdiff = 0.0;
    /* Map and partial reduce */
    while( loc < NUMLOCI){
        double elem = mapLogDiffsFunc(Q,TestQ,P,Geno,ind,loc);
        logdiff += elem;
        loc += get_global_size(0);
    }

    /* reduce locally */
    int localLoc = get_local_id(0);
    scratch[localLoc] = logdiff;
    barrier(CLK_LOCAL_MEM_FENCE);
    int devs = get_local_size(0);
    for(int offset = get_local_size(0) /2; offset > 0; offset >>= 1){
        if(localLoc < offset){
            scratch[localLoc] += scratch[localLoc + offset];
        }
        //Handle if were not working on a multiple of 2
        if (localLoc == 0 && (devs-1)/2 == offset){
            scratch[localLoc] += scratch[devs-1];
        }
        devs >>= 1;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    /* save result */
    int gid = get_group_id(0);
    if(localLoc == 0){
        results[ind*numgroups +gid] = scratch[0];
    }

    /* reduce over the groups into final result */
    barrier(CLK_GLOBAL_MEM_FENCE);
    if(gid==0){
        logdiffs[ind] = 0;
        for(int id =0; id < numgroups; id ++){
            logdiffs[ind] += results[ind*numgroups + id];
            results[ind*numgroups + id] = 0;
        }
    }
}

/* copy of logdiffs, but calc only one item */
double mapLogLikeFunc(__global double *Q, __global double *P,
                       __global int *Geno, int ind, int loc)
{
    int allele, line, pop;
    double term;
    if (ind < NUMINDS && loc < NUMLOCI){
        double runningtotal = 1.0;
        double logterm = 0.0;
        for (line = 0; line < LINES; line++) {
            allele = Geno[GenPos (ind, line, loc)];
            if (allele != MISSING) {
                term = 0.0;
                for (pop = 0; pop < MAXPOPS; pop++) {
                    term += Q[QPos(ind,pop)] * P[PPos (loc, pop, allele)];
                }

                //TODO: Evaluate underflow safe vs nonsafe
                // safe version, should not underflow
                
                if (term > SQUNDERFLO) {
                    runningtotal *= term;
                } else {
                    runningtotal *= SQUNDERFLO;
                }
                if (runningtotal < SQUNDERFLO){
                    logterm += log(runningtotal);
                    runningtotal = 1.0;
                }
                //Might underflow?
                /*logterm += log(term);*/
            }
        }
        logterm += log(runningtotal);
        return logterm;
    }
    return 0.0;
}

__kernel void mapReduceLogLike(__global double *Q,
                                __global double *P,
                                __global int *Geno,
                                __global double *loglikes,
                                __global double *results,
                                __local  double *scratch)
{
    int ind = get_global_id(1);
    int loc = get_global_id(0);
    int numgroups = get_num_groups(0);
    /* idempotent */
    double logterm = 0.0;
    /* Map and partial reduce */
    /* clear results buffer */
    for(int id =0; id < numgroups; id ++){
        results[ind*numgroups + id] = 0;
    }
    barrier(CLK_GLOBAL_MEM_FENCE);
    while( loc < NUMLOCI){
        double elem = mapLogLikeFunc(Q,P,Geno,ind,loc);
        logterm += elem;
        loc += get_global_size(0);
    }

    /* reduce locally */
    int localLoc = get_local_id(0);
    scratch[localLoc] = logterm;
    barrier(CLK_LOCAL_MEM_FENCE);
    int devs = get_local_size(0);
    for(int offset = get_local_size(0) /2; offset > 0; offset >>= 1){
        if(localLoc < offset){
            scratch[localLoc] += scratch[localLoc + offset];
        }
        //Handle if were not working on a multiple of 2
        if (localLoc == 0 && (devs-1)/2 == offset){
            scratch[localLoc] += scratch[devs-1];
        }
        devs >>= 1;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    /* save result */
    int gid = get_group_id(0);
    if(localLoc == 0){
        /*results[ind*numgroups +gid] = 1;*/
        results[ind*numgroups +gid] = scratch[0];
    }

    /* reduce over the groups into final result */
    barrier(CLK_GLOBAL_MEM_FENCE);
    if(gid==0){
        loglikes[ind] = 0;
        for(int id =0; id < numgroups; id ++){
            loglikes[ind] += results[ind*numgroups + id];
            results[ind*numgroups + id] = 0;
        }
    }
}

__kernel void CalcLike(
        __global double *loglikes,
        __global double *indlike_norm,
        __global double *sumindlike,
        const int usesumindlike, /* this is true after burnin */
        __global double *loglike,
        __global double *results,
        __local  double *scratch
        )
{
    int ind = get_global_id(0);

    double logterm = 0.0;
    int numgroups = get_num_groups(0);

    /* Map and partial reduce */
    while( ind < NUMINDS){
        double elem = loglikes[ind];
        if (usesumindlike) {
            if (indlike_norm[ind]==0.0) {
                indlike_norm[ind] = elem;
            }
            sumindlike[ind] += exp(elem-indlike_norm[ind]);
        }
        logterm += elem;
        ind += get_global_size(0);
    }

    int localLoc = get_local_id(0);
    scratch[localLoc] = logterm;
    barrier(CLK_LOCAL_MEM_FENCE);
    int devs = get_local_size(0);
    for(int offset = get_local_size(0) /2; offset > 0; offset >>= 1){
        if(localLoc < offset){
            scratch[localLoc] += scratch[localLoc + offset];
        }
        //Handle if were not working on a multiple of 2
        if (localLoc == 0 && (devs-1)/2 == offset){
            scratch[localLoc] += scratch[devs-1];
        }
        devs >>= 1;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    /* save result */
    int gid = get_group_id(0);
    if(localLoc == 0){
        results[gid] = scratch[0];
    }

    /* reduce over the groups into final result */
    barrier(CLK_GLOBAL_MEM_FENCE);
    if(gid==0){
        loglike[0] = 0;
        for(int id =0; id < numgroups; id++){
            loglike[0] += results[id];
            results[id] = 0;
        }
    }
}
