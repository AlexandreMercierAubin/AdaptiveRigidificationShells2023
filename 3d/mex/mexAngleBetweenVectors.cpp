/*==========================================================
 * computes the angle between two vectors wrt a rotation axis
 *
 * [ angles ] = mexAngleBetweenVectors( N1, N2, RotationAxis );
 *
 * inputs:
 *      N1           is n by 3 normals
 *      N2           is n by 3 normals
 *      RotationAxis is n by 3 vector
 *
 * outputs:
 *      angles       angles

 * To compile type: mex -R2018a mexAngleBetweenVectors.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    if ( nrhs != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:nrhs","three inputs required.");
    }
    // make sure the arguments are dense real double
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    if ( mxGetN(prhs[0]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:mustHave3cols","Parameters must have 3 columns.");
    }
    if ( mxGetN(prhs[1]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:mustHave3cols","Parameters must have 3 columns.");
    }
    if ( mxGetN(prhs[2]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAngleBetweenVectors:mustHave3cols","Parameters must have 3 columns.");
    }
    
    double *N1 = mxGetDoubles(prhs[0]);
    double *N2 = mxGetDoubles(prhs[1]);
    double *Axis = mxGetDoubles(prhs[2]);
    size_t n = mxGetM(prhs[0]); // could check that the n's match
    plhs[0] = mxCreateDoubleMatrix( n, 1, mxREAL );
    double *theta = mxGetDoubles(plhs[0]);
    
    for(size_t i = 0; i < n; ++i) {
        double N11 = N1[i];
        double N12 = N1[i+n];
        double N13 = N1[i+2*n];
        double N21 = N2[i];
        double N22 = N2[i+n];
        double N23 = N2[i+2*n];
        double A1 = Axis[i];
        double A2 = Axis[i+n];
        double A3 = Axis[i+2*n];
        double nc1 = A2*N13 - A3*N12;
        double nc2 = A3*N11 - A1*N13;
        double nc3 = A1*N12 - A2*N11;
        double ncm = sqrt( nc1*nc1 + nc2*nc2 + nc3*nc3 );
        double c1 = N11*N21 + N12*N22 + N13*N23;
        double c2 = (N21*nc1 + N22*nc2 + N23*nc3)/ ncm;        
		theta[i] = atan2(c2,c1);

// 		const Eigen::Vector3d n1 = N1.row(i);
// 		const Eigen::Vector3d n2 = N2.row(i);
// 		const Eigen::Vector3d r = RotationAxis.row(i);
// 		
// 		const Eigen::Vector3d ncross = r.cross(n1);
// 		const Eigen::Vector3d ncrossnorm = ncross/ncross.norm();
// 		
// 		Eigen::Matrix3d E;
// 		E.row(0) = r;
// 		E.row(1) = n1;
// 		E.row(2) = ncrossnorm;
// 		const Eigen::Vector3d c = E * n2;
// 		theta(i,0) = atan2(c[2],c[1]);
    }	
}
