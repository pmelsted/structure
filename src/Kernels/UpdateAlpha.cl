double AlphaPriorDiff (double newalpha, double oldalpha)
{
    /*returns log diff in priors for the alpha, assuming a gamma prior on alpha
      See notes 7/29/99 */
    return ((ALPHAPRIORA - 1) * log (newalpha / oldalpha) +
            (oldalpha - newalpha) / ALPHAPRIORB);
}


__kernel void UpdateAlpha(
       __global double *Q,
       __global double *Alpha,
       __global int *popflags,
       __global double *norms,
       __global double *results,
       __global uint *randGens,
       __local double *scratch,
       const int POPFLAGINDS)
{
}
/*{
    int ind = get_global_id(0);
    int pop = get_global_id(1);
    int numgroups = get_num_groups(0);
    int gid = get_group_id(0);

    int redpop;
    int numredpops = MAXPOPS;

    double newalpha = norms[pop];
    double oldalpha = Alpha[pop];
    double alphasum =0.0;

    double sum = 0.0;

    if ((newalpha > 0) && ((newalpha < ALPHAMAX) || (!(UNIFPRIORALPHA)) ) ) {
        for (int i=0; i<MAXPOPS; i++)  {
            alphasum += Alpha[i];
        }
        if (POPALPHAS) numredpops = pop +1;
        while( ind < NUMINDS){
            if (!((USEPOPINFO) && (popflags[ind]))) {
                double elem = 0.0;
                for(redpop = pop; redpop < numredpops; redpop++){
                    elem += log(Q[QPos (ind, redpop)]);
                }
                sum += elem;
                ind += get_global_size(0);
            }
        }

        [> reduce locally <]
        int localId = get_local_id(0);
        scratch[localId] = sum;
        barrier(CLK_LOCAL_MEM_FENCE);
        int devs = get_local_size(0);
        for(int offset = get_local_size(0) /2; offset > 0; offset >>= 1){
            if(localId < offset){
                scratch[localId] += scratch[localId + offset];
            }
            //Handle if were not working on a multiple of 2
            if (localId == 0 && (devs-1)/2 == offset){
                scratch[localId] += scratch[devs-1];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        if(localId == 0){
            results[pop*numgroups +gid] = scratch[0];
        }
        //TODO: Handle if numgroups are more than MAXGROUPS
        //Possibly by reducing with global barrier.

        barrier(CLK_GLOBAL_MEM_FENCE);
        if(gid==0){
            RndDiscState randState[1];
            initRndDiscState(randState,randGens,pop);
            [>int multiple = 1;<]
            [>if (ONEFST) multiple = MAXPOPS;<]
            [>double logprobdiff = FPriorDiff (newf, oldf);<]
            [>logprobdiff += multiple*NUMLOCI*lgamma(newfrac);<]
            [>logprobdiff -= multiple*NUMLOCI*lgamma(oldfrac);<]
            double logprobdiff = 0.0;
            double logterm = 0.0;
            if (!(UNIFPRIORALPHA)) logprobdiff = AlphaPriorDiff (newalpha, oldalpha);

            for(int id =0; id < numgroups; id ++){
                logterm += results[pop*numgroups + id];
            }

            int multiple = numredpops - pop;
            double lpsum = 0.0;
            lpsum -= (oldalpha - 1.0) * logterm;
            lpsum += (newalpha - 1.0) * logterm;

            double sumalphas = alphasum;
            lpsum -= (lgamma (sumalphas) - multiple * lgamma ( oldalpha)) * POPFLAGINDS;

            if (POPALPHAS){
                sumalphas += newalpha - oldalpha;
            } else {
                sumalphas = MAXPOPS*newalpha;
            }

            lpsum += (lgamma (sumalphas) - multiple * lgamma ( newalpha)) * POPFLAGINDS;
            logprobdiff += lpsum;

            if (rndDisc(randState) < exp(logprobdiff)) {   [>accept new f <]
                for(redpop = pop; redpop < numredpops; redpop++){
                    Alpha[redpop] = newalpha;
                }
            }
            saveRndDiscState(randState);
        }
    }
}*/
