__kernel void NonIndUpdateEpsilon(
        __global float *P,
        __global float *Epsilon,
        __global float *Fst,
        __global int *NumAlleles,
        __global uint *randGens,
        __global float *lambdas,
        const float invsqrtnuminds)
{

    int loc = get_global_id(0);
    int allele1,allele2;
    int eps1,eps2;
    int pop;
    float diff;
    float sum;
    float lambda = lambdas[0];
    float rand;
    RndDiscState randState[1];

    while (loc < NUMLOCI){
        if (NumAlleles[loc] > 1){
            initRndDiscState(randState,randGens,loc);

            allele1 = RandomInteger(0,NumAlleles[loc],randState);
            allele2 = RandomInteger(0,NumAlleles[loc]-1,randState);
            if (allele2 >= allele1){
                allele2 += 1;
            }
            eps1 = Epsilon[EpsPos(loc,allele1)];
            eps2 = Epsilon[EpsPos(loc,allele2)];
            diff = numToRange(0,invsqrtnuminds,rndDisc(randState));

            if(eps1 + diff < 1.0 && eps2 - diff > 0.0){
                //TODO: Evaluate whether we should reduce here.
                sum=0.0;
                for (pop=0; pop<MAXPOPS; pop++) { /*compute likelihood ratio*/
                    float frac = (1.0-Fst[pop])/Fst[pop];
                    sum += lgamma(frac*eps1);
                    sum += lgamma(frac*eps2);
                    sum -= lgamma(frac*(eps1+diff));
                    sum -= lgamma(frac*(eps2-diff));

                    sum += frac*diff*log(P[PPos(loc, pop, allele1)]);
                    sum -= frac*diff*log(P[PPos(loc, pop, allele2)]);
                }
                if (lambda != 1.0) {              /*compute prior ratio*/
                    /* as it is in code */
                    /*float ratio = (eps1 + diff)* (eps2 - diff)/(eps1)/(eps2)*/
                    /*sum += log(pow(ratio, lambda-1.0));*/
                    /* as it probably should be ? */
                    float ratio = (eps1 + diff)* (eps2 - diff)/(eps1*eps2);
                    sum += (lambda-1.0)*log(ratio);
                }

                if (rndDisc(randState) < exp(sum)){
                    Epsilon[EpsPos(loc,allele1)]+=diff;
                    Epsilon[EpsPos(loc,allele2)]-=diff;
                }
            }
            saveRndDiscState(randState);
        }
        loc += get_global_size(0);
        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}
