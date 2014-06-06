#include "mwc64x.cl"
#define RAND_MAX 4294967296.0f

/*
 * returns a random real in the [lower,upper)
 */
float randomReal(float lower, float upper,mwc64x_state_t *randState){
    uint randVal;
    float randPercent;
    float random;

    randVal = MWC64X_NextUint(randState);
    randPercent = (float) randVal/(RAND_MAX +1);
    return (lower + randPercent*(upper-lower));
}

double numToRange(double low, double high, double num)
{
  /* Takes a number in [0,1) -> [low,high) */
  return (low + num * (high - low) );
}

int dimLoc(int * dims, int * dimMaxs, int numDims){
    int loc = 0;
    int i, j;
    for(i = 0; i < numDims;++i){
        int dimProd =1;
        for(j = i+1; j < numDims; j++){
            dimProd *= dimMaxs[j];
        }
        loc += dims[i]*dimProd;
    }
    return loc;
}
/*
 * Copies the entire last dimension over into localarr
 */
void copyToLocal( __global double * globalArr, double *localArr,
                  int * dims, int * dimMaxs, int numDims){
    int i;
    int origLastDim = dims[numDims-1];
    for(i = 0; i < dimMaxs[numDims-1]; ++i){
        dims[numDims-1] = i;
        localArr[i] = globalArr[dimLoc(dims,dimMaxs,numDims)];
    }
    dims[numDims-1] = origLastDim;
}

double rnd(double * localRandom, int * randomValsTaken){
    return localRandom[(*randomValsTaken)++];
    /*double value;*/
    /*do {*/
        /*value = localRandom[*randomValsTaken];*/
        /*(*randomValsTaken)++;*/
    /*} while ((value == 0.0 || value == 1.0));*/
    /*return value;*/
}


/*
 *  Returns a random number between 0 and n-1, according to a list of
 *  probabilities.  The function takes a (possibly) unnormalized
 *  vector of probabilities of the form (p1, p2, p3,....,p_n).  There
 *  are "total" possible options, and sum is the sum of the
 *  probabilities.  This comes up in the Gibbs sampler context.
 */
int PickAnOptionDiscrete(int total, double sum, double Probs [], double randNum){
   int option;
   double random;
   double sumsofar =  0.0;

   random = numToRange(0,sum, randNum);
   for (option=0; option < total; ++option){
       sumsofar += Probs[option];
	   if (random <= sumsofar) break;
   }
   return option;
}
