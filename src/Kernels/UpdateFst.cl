double
FPriorDiff (double newf, double oldf)
{
    /*returns log diff in priors for the correlation, f. See notes 5/14/99, and 7/15/99 */

    return ((FPRIORMEAN*FPRIORMEAN/(FPRIORSD*FPRIORSD) - 1) * log (newf / oldf) +
            (oldf - newf) *FPRIORMEAN/(FPRIORSD*FPRIORSD));
}






__kernel void UpdateFST(
            __global double *Epsilon,
            __global double *Fst,
            __global double *P,
            __global int *NumAlleles,
            __global double *normals,
            __global uint *randGens
            )
{
    int loc = get_global_id(0);
    int pop = get_global_id(1);

    double oldf = Fst[pop];
    double newf = normals[pop];
    if (newf > 0.0 && newf < 1.0){

    }

}
