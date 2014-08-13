float AlphaPriorDiff (float newalpha, float oldalpha)
{
    /*returns log diff in priors for the alpha, assuming a gamma prior on alpha
      See notes 7/29/99 */
    return ((ALPHAPRIORA - 1) * log (newalpha / oldalpha) +
            (oldalpha - newalpha) / ALPHAPRIORB);
}


__kernel void UpdateAlpha(
       __global float *Q,
       __global float *Alpha,
       __global int *popflags,
       __global float *norms,
       __global float *results,
       __global uint *randGens,
       __local float *scratch,
       const int POPFLAGINDS)
{
    int alpha = get_global_id(1);
    while( alpha < NUMALPHAS){
        int ind = get_global_id(0);
        if(ind < NUMINDS){
            int redpop;
            int numredpops = MAXPOPS;

            float newalpha = norms[alpha];
            float oldalpha = Alpha[alpha];
            float alphasum =0.0;

            if ((newalpha > 0) && ((newalpha < ALPHAMAX) || (!(UNIFPRIORALPHA)) ) ) {
                if (POPALPHAS){ numredpops = alpha +1; }
                //TODO: Evaluate underflow safe vs nonsafe
                float sum = 1.0;
                float total = 0.0;
                while( ind < NUMINDS){
                    if (!((USEPOPINFO) && (popflags[ind]))) {
                        //Safe version (similar to in code)
                        //Watching out for underflow
                        float elem = 1.0;
                        for(redpop = alpha; redpop < numredpops; redpop++){
                            elem *= Q[QPos (ind, redpop)];
                        }

                        if (elem > SQUNDERFLO){
                            sum *= elem;
                        } else {
                            sum *= SQUNDERFLO;
                        }
                        if(sum < SQUNDERFLO){
                            total += log(sum);
                            sum = 1.0;
                        }
                        //Might underflow?
                      /*float elem = 0.0;
                        for(redpop = alpha; redpop < numredpops; redpop++){
                            elem += log(Q[QPos (ind, redpop)]);
                        }

                        total += elem;
                      */
                        ind += get_global_size(0);
                    }
                }

                total += log(sum);

                /* reduce locally */
                int localId = get_local_id(0);
                scratch[localId] = total;
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

                    devs >>= 1;
                    barrier(CLK_LOCAL_MEM_FENCE);
                }

                int numgroups = get_num_groups(0);
                int gid = get_group_id(0);

                if(localId == 0){
                    results[alpha*numgroups +gid] = scratch[0];
                }
                //TODO: Handle if numgroups are more than MAXGROUPS
                //Possibly by reducing with global barrier.

                barrier(CLK_GLOBAL_MEM_FENCE);
                if(gid==0){
                    RndDiscState randState[1];
                    initRndDiscState(randState,randGens,alpha);
                    for (int i=0; i<MAXPOPS; i++)  {
                        alphasum += Alpha[i];
                    }
                    float logprobdiff = 0.0;
                    float logterm = 0.0;
                    if (!(UNIFPRIORALPHA)) logprobdiff = AlphaPriorDiff (newalpha, oldalpha);

                    for(int id =0; id < numgroups; id++){
                        logterm += results[alpha*numgroups + id];
                        results[alpha*numgroups + id] = 0;
                    }

                    int multiple = numredpops - alpha;
                    float lpsum = (newalpha - oldalpha) * logterm;
                    /*lpsum -= (oldalpha - 1.0) * logterm;
                      lpsum += (newalpha - 1.0) * logterm;*/

                    float sumalphas = alphasum;
                    if (POPALPHAS){
                        sumalphas += newalpha - oldalpha;
                    } else {
                        sumalphas = MAXPOPS*newalpha;
                    }

                    lpsum += ((lgamma(sumalphas) - lgamma(alphasum)) - multiple*(lgamma(newalpha) - lgamma(oldalpha))) * POPFLAGINDS;
                    
                    /*lpsum -= (lgamma (alphasum) - multiple * lgamma ( oldalpha)) * POPFLAGINDS;
                      lpsum += (lgamma (sumalphas) - multiple * lgamma ( newalpha)) * POPFLAGINDS;*/
                    logprobdiff += lpsum;

                    if (rndDisc(randState) < exp(logprobdiff)) {   /*accept new f */
                        for(redpop = alpha; redpop < numredpops; redpop++){
                            Alpha[redpop] = newalpha;
                        }
                    }
                    saveRndDiscState(randState);
                }
            }
        }
        alpha += get_global_size(1);
        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}
