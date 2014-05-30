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


/*
 *  Returns a random number between 0 and n-1, according to a list of
 *  probabilities.  The function takes a (possibly) unnormalized
 *  vector of probabilities of the form (p1, p2, p3,....,p_n).  There
 *  are "total" possible options, and sum is the sum of the
 *  probabilities.  This comes up in the Gibbs sampler context.
 */
int PickAnOption(int total, float sum, float Probs [], mwc64x_state_t *randState){
   int option;
   float random;
   float sumsofar;
   
   sumsofar = 0.0; 
   random = randomReal(0,sum, randState);
   for (option=0; option < total && random <= sumsofar; ++option){
       sumsofar += Probs[option];
   }
   return option;
}
