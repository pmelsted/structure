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
       __global uint *randGens)
{
    int loc = get_global_id(0);
    int pop = get_global_id(1);
    int numgroups = get_num_groups(0);
    int gid = get_group_id(0);

    int redpop;
    int numredpops;

    double newalpha = norms[pop];
    double oldalpha = Alpha[pop];

    if ((newalpha > 0) && ((newalpha < ALPHAMAX) || (!(UNIFPRIORALPHA)) ) ) {

        numredpops = pop +1;
        if (!(POPALPHAS)) numredpops = MAXPOPS;
        //TODO: Handle if numgroups are more than MAXGROUPS
        //Possibly by reducing with global barrier.

        barrier(CLK_GLOBAL_MEM_FENCE);
        if(gid==0){
            RndDiscState randState[1];
            initRndDiscState(randState,randGens,pop);
            /*int multiple = 1;*/
            /*if (ONEFST) multiple = MAXPOPS;*/
            /*double logprobdiff = FPriorDiff (newf, oldf);*/
            /*logprobdiff += multiple*NUMLOCI*lgamma(newfrac);*/
            /*logprobdiff -= multiple*NUMLOCI*lgamma(oldfrac);*/
            double logprobdiff = 0.0;
            if (!(UNIFPRIORALPHA)) logprobdiff = AlphaPriorDiff (newalpha, oldalpha);

            for(int id =0; id < numgroups; id ++){
                logprobdiff += results[pop*numgroups + id];
            }
            if (rndDisc(randState) < exp(logprobdiff)) {   /*accept new f */
                for(redpop = pop; redpop < numredpops; redpop++){
                    Alpha[redpop] = newalpha;
                }
            }
            saveRndDiscState(randState);
        }
    }
}
