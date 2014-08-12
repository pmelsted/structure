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

float numToRange(float low, float high, float num)
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
void copyToLocal( __global float * globalArr, float *localArr,
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
int PickAnOptionDiscrete(int total, float sum, float Probs [],
                         RndDiscState *randState)
{
    int option;
    float random;
    float sumsofar =  0.0;

    random = numToRange(0,sum, rndDisc(randState));
    for (option=0; option < total; ++option) {
        sumsofar += Probs[option];
        if (random <= sumsofar) {
            break;
        }
    }
    return option;
}

/* int RandomInteger(int low, int high,RndDiscState *randState) */
/* { */
/*     if (high == low){ */
/*         return low; */
/*     } */
/*     int range = high-low; */
/*     uint random = rndUInt(randState) % range; */
/*     return (int) random + low; */
/* } */

int RandomInteger(int low, int high,RndDiscState *randState)
{
    int k;
    float d = rndDisc(randState);

    k = (int) (d * (high - low + 1));
    return (low + k);
}


int AlphaPos(int loc, int pop)
{
    if((loc)==NUMLOCATIONS) {
        return pop;
    } else {
        return (MAXPOPS*((loc)+1)+(pop));
    }
}

/* Returns gamma(f,1), where 0 < f < 1 */
float RGammaDiscFloat(float n,RndDiscState *randState){
    float x=0.0;
    float E=2.71828182;
    float b=(n+E)/E;
    float p=0.0;
    while(1){
        p=b*rndDisc(randState);
        if(p>1) {
            x=-log((b-p)/n);
            if (!(((n-1)*log(x))<log(rndDisc(randState)))) {
                /* Accept */
                break;
            }
        } else {
            x=exp(log(p)/n);
            if(!(x>-log(rndDisc(randState)))) {
                /* Accept */
                break;
            }
        }
    }
    return x;
}

/* Returns gamma(1,1) */
float RGammaDiscOne(RndDiscState *randState){
    float a=0.0;
    float u,u0,ustar;
    u=rndDisc(randState);
    u0=u;
    while (1){
        ustar=rndDisc(randState);
        if(u<ustar) {
            break;
        }
        u=rndDisc(randState);
        if(u>=ustar) {
            a += 1;
            u=rndDisc(randState);
            u0=u;
        }
    }
    return (a+u0);
}

/* Returns gamma(n,1) where n is an int */
float RGammaDiscInt(int n,RndDiscState *randState){
    int i =0;
    float x = 0;
    for(i = 0; i < n; ++i){
        x += log(rndDisc(randState));
    }
    return -x;
}

/*  taken from page 413 of
 *
 *  Non-Uniform Random Variate Generation by Luc Devroye
 *
 *  (originally published with Springer-Verlag, New York, 1986)
 *
 *  which can be found online at http://luc.devroye.org/rnbookindex.html
 */
float RGammaCheng(float a,RndDiscState *randState){
    float b = a - log(4.0);
    float c = a + sqrt(2.0*a-1.0);
    float U,V,X,Y,Z,R;
    while (1){
        U = rndDisc(randState);
        V = rndDisc(randState);
        Y = a*log(V/(1.0-V));
        X = a*exp(V);
        Z = U*(V*V);
        R = b + c*Y - X;
        if( (R >= (9.0/2.0)*Z - (1+log(9.0/2.0))) ||   ( R >= log(Z)) ){
            break;
        }
    }
    return X;
}

/*  taken from page 410 of
 *
 *  Non-Uniform Random Variate Generation by Luc Devroye
 *
 *  (originally published with Springer-Verlag, New York, 1986)
 *
 *  which can be found online at http://luc.devroye.org/rnbookindex.html
 */

float RGammaBest(float a,RndDiscState *randState){
    float b = a -1;
    float c = 3*a - 0.75;
    float U,V,W,Y,X,Z;
    while (1){
        U = rndDisc(randState);
        V = rndDisc(randState);
        W = U*(1-U);
        Y = sqrt((c/W))*(U-0.5);
        X = b + Y;
        if( X >= 0){
            Z = 64*(W*W*W)*(V*V);
            if((Z <= 1 - (2*(Y*Y))/X) || (log(Z) <= 2*(b*log(X/b)-Y))){
                break;
            }
        }
    }
    return X;
}

float RGammaLargerThanOne(float n, RndDiscState *randState)
{
    float aa,w,x;
    float nprev=0.0;
    float c1=0.0;
    float c2=0.0;
    float c3=0.0;
    float c4=0.0;
    float c5=0.0;
    float u1;
    float u2;
    /*if(n!=nprev) {*/
        c1=n-1.0;
        aa=1.0/c1;
        c2=aa*(n-1/(6*n));
        c3=2*aa;
        c4=c3+2;
        if(n>2.5) {
            c5=1/sqrt(n);
        }
    /*}*/
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
    return x;

}

/*-----------Gamma and dirichlet from Matt.----------*/
/* gamma random generator from Ripley, 1987, P230 */

float RGammaDisc(float n,float lambda,RndDiscState *randState)
{
    float x=0.0;
    if(n<1) {
        x = RGammaDiscFloat(n,randState);
    } else if(n==1.0) {
        /* gamma (1,1) is an exponential dist */
        /*x = RGammaDiscOne(randState);*/
        x = -log(rndDisc(randState));
    } else {
        /*if (rndDisc(randState) < 0.5){*/
        x = RGammaBest(n,randState);
        /* x = RGammaLargerThanOne(n,randState); */
        /* x = RGammaCheng(n,randState); */
        /*} else {*/
            /*x = RGammaCheng(n,randState);*/
        /*}*/
        /*int wholepart = (int) n;*/
        /*float xi = 0.0;*/
        /*float wholegamma = 0.0;*/
        /*float floatpart = n - wholepart;*/
        /*xi = RGammaDiscFloat(floatpart,randState);*/
        /*wholegamma = RGammaBest(wholepart,randState);*/
        /*x = xi + wholegamma;*/
    }
    return x/lambda;
}





/* Dirichlet random generator
   a and b are arrays of length k, containing floats.
   a is the array of parameters
   b is the output array, where b ~ Dirichlet(a)
   */

void RDirichletDisc(float * a, int k, float * b,RndDiscState *randState)
{
    int i;
    float sum=0.0;
    for(i=0; i<k; i++) {
        b[i]=RGammaDisc(a[i],1,randState);
        sum += b[i];
    }
    for(i=0; i<k; i++) {
        b[i] /= sum;
    }
}

/* Melissa's version, adapted from an algorithm on wikipedia.  January 08 */
float LogRGammaDisc(float n, float lambda, RndDiscState *randState)
{
    float v0, v[3], E=2.71828182, em, logem, lognm;
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
void LogRDirichletDisc (float *a, int k,__global float *b,
                        RndDiscState *randState)
{
    int i;
    float sum = 0.0;
    /* float sum2; */
    for (i = 0; i < k; i++) {
        b[i] =RGammaDisc (a[i], 1,randState);
        sum += b[i];
    }

    /* patch added May 2007 to set gene frequencies equal if all draws
       from the Gamma distribution are very low.
       Ensures that P and logP remain defined in this rare event */
    /* if(sum<UNDERFLO) { */
    /*     for(i=0; i<k; i++) { */
    /*         b[i] = 1.0/(float)(k); */
    /*     } */
    /* } else { */
        /* sum2=log(sum); */
        for (i = 0; i < k; i++) {
            b[i]/=sum;
        }
    /* } */
}
