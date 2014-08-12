double
FPriorDiff (double newf, double oldf)
{
    /*returns log diff in priors for the correlation, f. See notes 5/14/99, and 7/15/99 */
    double pri = (FPRIORMEAN*FPRIORMEAN/(FPRIORSD*FPRIORSD) - 1);

    return (pri * (log (newf) - log( oldf)) + (oldf - newf) *FPRIORMEAN/(FPRIORSD*FPRIORSD));
}



double FlikeFreqsDiffMap (double newfrac,double oldfrac,
        __global double *Epsilon,
        __global double *P,
        __global int *NumAlleles,
        int loc,int pop){
    int allele;
    double eps,logp;
    double sum;

    if (NumAlleles[loc]==0) {
        return -(lgamma(newfrac) - lgamma(oldfrac)); /* should not be counting sites with all missing data */
    } else {
        sum = 0.0;
        for (allele=0; allele < NumAlleles[loc]; allele++) {
            eps = Epsilon[EpsPos (loc, allele)];
            logp = log(P[PPos(loc,pop,allele)]);
            sum += newfrac*eps*logp;
            sum -= oldfrac*eps*logp;
            sum -= lgamma( newfrac*eps);
            sum += lgamma( oldfrac*eps);
        }
        return sum;
    }
}



__kernel void UpdateFst(
            __global double *Epsilon,
            __global double *Fst,
            __global double *P,
            __global int *NumAlleles,
            __global double *normals,
            __global uint *randGens,
            __global double *results,
            __local  double *scratch)
{
    int pop = get_global_id(1);
    int numgroups = get_num_groups(0);
    while (pop < MAXPOPS){
        int loc = get_global_id(0);
        double newf = normals[pop];
        /* ensure newf is large enough so we don't cause over/underflow */
        if (newf > 10e-10 && newf < 1.0){
            double sum = 0.0;
            int redpop;
            int numredpops;
            double oldf = Fst[pop];
            double newfrac = (1.0-newf)/newf;
            double oldfrac = (1.0-oldf)/oldf;
            numredpops = pop +1;
            if (ONEFST) numredpops = MAXPOPS;
            /* idempotent */
            /* Map and partial reduce */
            while( loc < NUMLOCI){
                double elem = 0.0;
                for(redpop = pop; redpop < numredpops; redpop++){
                    elem += FlikeFreqsDiffMap(newfrac,oldfrac,Epsilon,P,NumAlleles,loc,redpop);
                }
                sum += elem;
                loc += get_global_size(0);
            }
            /* reduce locally */
            int localLoc = get_local_id(0);
            scratch[localLoc] = sum;
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

            int gid = get_group_id(0);
            //TODO: Handle if numgroups are more than MAXGROUPS
            //Possibly by reducing with global barrier.
            if(localLoc == 0){
                results[pop*numgroups +gid] = scratch[0];
            }

            barrier(CLK_GLOBAL_MEM_FENCE);
            if(gid==0){
                RndDiscState randState[1];
                initRndDiscState(randState,randGens,pop);
                int multiple = 1;
                if (ONEFST) multiple = MAXPOPS;
                double logprobdiff = FPriorDiff (newf, oldf);
                logprobdiff += multiple*NUMLOCI*lgamma(newfrac);
                logprobdiff -= multiple*NUMLOCI*lgamma(oldfrac);
                for(int id =0; id < numgroups; id ++){
                    logprobdiff += results[pop*numgroups + id];
                    results[pop*numgroups + id] = 0;
                }
                if (logprobdiff >= 0.0 || rndDisc(randState) < exp(logprobdiff)) {   /*accept new f */
                    for(redpop = pop; redpop < numredpops; redpop++){
                        Fst[redpop] = newf;
                    }
                }
                saveRndDiscState(randState);
            }
        }
        pop += get_global_size(1);
        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}

