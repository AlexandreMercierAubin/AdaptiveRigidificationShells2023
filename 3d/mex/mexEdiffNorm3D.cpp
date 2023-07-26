/*==========================================================
 * mexEdiffNorm.cpp - computes 3D strain difference Frobeneneous norm squared
 *
 * [ eDotNormsSquared ] = mexEdotNorm( Fa, Fb, h );
 *
 * Input:
 *  Fa      9x#E deformation gradient of each element (or a vector of size 9x#E) at one time step
 *  Fb      9x#E deformation gradient of each element (or a vector of size 9x#E) at another time step
 *  h       timestep
 *
 * Output:
 *  EdiffNormSquared     #E Frobeneus norm squared of strain difference 
 *                          divided by h), i.e., an approximation of strain rate
 *
 * To compile type: mex -R2018a mexEdiffNorm3D.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *Fa = mxGetDoubles(prhs[0]);
    size_t sizeF = mxGetM(prhs[0]);
    double *Fb = mxGetDoubles(prhs[1]);
    double h = mxGetScalar( prhs[2] );

    if ( nrhs != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexEdiffNorm:nrhs","3 inputs required.");
    }
    // make sure the arguments are dense real double
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ) {
        mexErrMsgIdAndTxt("ARP:mexEdiffNorm:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) ) {
        mexErrMsgIdAndTxt("ARP:mexEdiffNorm:mustBeDenseRealDouble","Parameter must be dense real double");
    }

    if ( mxGetM(prhs[0]) != mxGetM(prhs[1]) || mxGetN(prhs[0]) != mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("ARP:mexEdiffNorm:equalSize","size of Fa should be equal to Fb");
    }

    if ( mxGetM(prhs[0])%9 != 0) {
        mexErrMsgIdAndTxt("ARP:mexEdiffNorm:power9","Fa should be a power of 9 entries");
    }

    if ( mxGetN(prhs[2]) != 1 ||  mxGetM(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("ARP:mexEdiffNorm:mustHave3cols","h should only be 1 by 1");
    }

    int sizeEdiffNorm = sizeF/9.0;
    plhs[0] = mxCreateDoubleMatrix( sizeEdiffNorm, 1, mxREAL );
    double *EdiffNormSquared = mxGetDoubles(plhs[0]);
    
    for(size_t i = 0; i < sizeEdiffNorm; ++i) {
        int pos = i*9;
        double Fa1_1 = F[pos];
        double Fa2_1 = F[pos+1];
        double Fa3_1 = F[pos+2];
        double Fa1_2 = F[pos+3];
        double Fa2_2 = F[pos+4];
        double Fa3_2 = F[pos+5];
        double Fa1_3 = F[pos+6];
        double Fa2_3 = F[pos+7];
        double Fa3_3 = F[pos+8];

        double Fb1_1 = Fd[pos];
        double Fb2_1 = Fd[pos+1];
        double Fb3_1 = Fd[pos+2];
        double Fb1_2 = Fd[pos+3];
        double Fb2_2 = Fd[pos+4];
        double Fb3_2 = Fd[pos+5];
        double Fb1_3 = Fd[pos+6];
        double Fb2_3 = Fd[pos+7];
        double Fb3_3 = Fd[pos+8];
        
        double t2 = 1.0/(h*h);
        double t0 = t2*pow((Fa1_1*Fa1_2)/2.0+(Fa2_1*Fa2_2)/2.0+(Fa3_1*Fa3_2)/2.0-(Fb1_1*Fb1_2)/2.0-(Fb2_1*Fb2_2)/2.0-(Fb3_1*Fb3_2)/2.0,2.0)*2.0+t2*pow((Fa1_1*Fa1_3)/2.0+(Fa2_1*Fa2_3)/2.0+(Fa3_1*Fa3_3)/2.0-(Fb1_1*Fb1_3)/2.0-(Fb2_1*Fb2_3)/2.0-(Fb3_1*Fb3_3)/2.0,2.0)*2.0+t2*pow((Fa1_2*Fa1_3)/2.0+(Fa2_2*Fa2_3)/2.0+(Fa3_2*Fa3_3)/2.0-(Fb1_2*Fb1_3)/2.0-(Fb2_2*Fb2_3)/2.0-(Fb3_2*Fb3_3)/2.0,2.0)*2.0+t2*pow((Fa1_1*Fa1_1)/2.0+(Fa2_1*Fa2_1)/2.0+(Fa3_1*Fa3_1)/2.0-(Fb1_1*Fb1_1)/2.0-(Fb2_1*Fb2_1)/2.0-(Fb3_1*Fb3_1)/2.0,2.0)+t2*pow((Fa1_2*Fa1_2)/2.0+(Fa2_2*Fa2_2)/2.0+(Fa3_2*Fa3_2)/2.0-(Fb1_2*Fb1_2)/2.0-(Fb2_2*Fb2_2)/2.0-(Fb3_2*Fb3_2)/2.0,2.0)+t2*pow((Fa1_3*Fa1_3)/2.0+(Fa2_3*Fa2_3)/2.0+(Fa3_3*Fa3_3)/2.0-(Fb1_3*Fb1_3)/2.0-(Fb2_3*Fb2_3)/2.0-(Fb3_3*Fb3_3)/2.0,2.0);

        EdiffNormSquared[i] = t0;
    }
}
