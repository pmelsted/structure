#include "mwc64x.cl"
#include "randGen.cl"
#define RAND_MAX 4294967296.0f

/*
 * returns a random real in the [lower,upper)
 */
float randomReal(float lower, float upper,mwc64x_state_t *randState)
{
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

int dimLoc(int * dims, int * dimMaxs, int numDims)
{
    int loc = 0;
    int i, j;
    for(i = 0; i < numDims; ++i) {
        int dimProd =1;
        for(j = i+1; j < numDims; j++) {
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
                  int * dims, int * dimMaxs, int numDims)
{
    int i;
    int origLastDim = dims[numDims-1];
    for(i = 0; i < dimMaxs[numDims-1]; ++i) {
        dims[numDims-1] = i;
        localArr[i] = globalArr[dimLoc(dims,dimMaxs,numDims)];
    }
    dims[numDims-1] = origLastDim;
}


/*
 *  Returns a random number between 0 and n-1, according to a list of
 *  probabilities.  The function takes a (possibly) unnormalized
 *  vector of probabilities of the form (p1, p2, p3,....,p_n).  There
 *  are "total" possible options, and sum is the sum of the
 *  probabilities.  This comes up in the Gibbs sampler context.
 */
int PickAnOptionDiscrete(int total, double sum, double Probs [],
                         RndDiscState *randState)
{
    int option;
    double random;
    double sumsofar =  0.0;

    random = numToRange(0,sum, rndDisc(randState));
    for (option=0; option < total; ++option) {
        sumsofar += Probs[option];
        if (random <= sumsofar) {
            break;
        }
    }
    return option;
}


int AlphaPos(int loc, int pop)
{
    if((loc)==NUMLOCATIONS) {
        return pop;
    } else {
        return (MAXPOPS*((loc)+1)+(pop));
    }
}

/*-----------Gamma and dirichlet from Matt.----------*/
/* gamma random generator from Ripley, 1987, P230 */

double RGammaDisc(double n,double lambda,RndDiscState *randState)
{

    double aa;
    double w;

    double x=0.0;
    if(n<1) {
        double E=2.71828182;
        double b=(n+E)/E;
        double p=0.0;
one:
        p=b*rndDisc(randState);
        if(p>1) {
            goto two;
        }
        x=exp(log(p)/n);
        if(x>-log(rndDisc(randState))) {
            goto one;
        }
        goto three;
two:
        x=-log((b-p)/n);
        if (((n-1)*log(x))<log(rndDisc(randState))) {
            goto one;
        }
three:
        ;
    } else if(n==1.0) {

        double a=0.0;
        double u,u0,ustar;
ten:
        u=rndDisc(randState);
        u0=u;
twenty:
        ustar=rndDisc(randState);
        if(u<ustar) {
            goto thirty;
        }
        u=rndDisc(randState);
        if(u<ustar) {
            goto twenty;
        }
        a += 1;
        goto ten;
thirty:
        return (a+u0)/lambda;
    } else {
        double nprev=0.0;
        double c1=0.0;
        double c2=0.0;
        double c3=0.0;
        double c4=0.0;
        double c5=0.0;
        double u1;
        double u2;
        if(n!=nprev) {
            c1=n-1.0;
            aa=1.0/c1;
            c2=aa*(n-1/(6*n));
            c3=2*aa;
            c4=c3+2;
            if(n>2.5) {
                c5=1/sqrt(n);
            }
        }
four:
        u1=rndDisc(randState);
        u2=rndDisc(randState);
        if(n<=2.5) {
            goto five;
        }
        u1=u2+c5*(1-1.86*u1);
        if ((u1<=0) || (u1>=1)) {
            goto four;
        }
five:
        w=c2*u2/u1;
        if(c3*u1+w+1.0/w < c4) {
            goto six;
        }
        if(c3*log(u1)-log(w)+w >=1) {
            goto four;
        }
six:
        x=c1*w;
        nprev=n;
    }

    return x/lambda;
}


/* Dirichlet random generator
   a and b are arrays of length k, containing doubles.
   a is the array of parameters
   b is the output array, where b ~ Dirichlet(a)
   */

void RDirichletDisc(double * a, int k, double * b,RndDiscState *randState)
{
    int i;
    double sum=0.0;
    for(i=0; i<k; i++) {
        b[i]=RGammaDisc(a[i],1,randState);
        sum += b[i];
    }
    for(i=0; i<k; i++) {
        b[i] /= sum;
    }
}

/* Melissa's version, adapted from an algorithm on wikipedia.  January 08 */
double LogRGammaDisc(double n, double lambda, RndDiscState *randState)
{
    double v0, v[3], E=2.71828182, em, logem, lognm;
    int i;
    if (n >= 1.0) {
        return log(RGammaDisc(n, lambda,randState));
    }
    v0 = E/(E+n);
    while (1) {
        for (i=0; i<3; i++) {
            v[i] = rndDisc(randState);
        }

        if (v[0] <= v0) {
            logem = 1.0/n*log(v[1]);
            em = exp(logem);
            lognm = log(v[2])+(n-1)*logem;
        } else {
            em = 1.0-log(v[1]);
            logem = log(em);
            lognm = log(v[2]) - em;
        }
        if (lognm <= (n-1)*logem - em) {
            return logem - log(lambda);
        }
    }
}

/*O(k)*/
void LogRDirichletDisc (double *a, int k,__global double *b,
                        RndDiscState *randState)
{
    int i;
    double sum = 0.0;
    double sum2;
    for (i = 0; i < k; i++) {
        b[i] =RGammaDisc (a[i], 1,randState);
        sum += b[i];
    }

    /* patch added May 2007 to set gene frequencies equal if all draws
       from the Gamma distribution are very low.
       Ensures that P and logP remain defined in this rare event */
    /*if(sum<UNDERFLO) {
        for(i=0; i<k; i++) {
            b[i] = 1.0/(double)(k);
            c[i] = log(b[i]);
        }
    } else {*/
    sum2=log(sum);
    for (i = 0; i < k; i++) {
        b[i]/=sum;
    }
    /*}*/
}
