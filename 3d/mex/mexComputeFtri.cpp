/* 
 * [ N ] = mexComputeFtri( V, T , DmInv);
 *
 * Input:
 *  V        #vertices by 3 matrix of vertex position
 *  T        #E by 3 vertices of elements
 * 
 * Output:
 *  F         #F 4el vector of 2 by 2 deformation gradients
 *
 * To compile type: mex -R2018a mexComputeFtri.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h> 

/* 
 * The gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    const size_t numDOFs = mxGetM(prhs[0]);
    double *V = mxGetDoubles(prhs[0]);
    const size_t numVerts = mxGetM(prhs[0]);
    int* T = (int*)mxGetData(prhs[1]);       //#tri by 3 indices (starting from 1) of verts of each element
    const size_t numElements = mxGetM(prhs[1]);
    double* DmInv = mxGetDoubles(prhs[2]);

    // create the output vectors 
    plhs[0] = mxCreateDoubleMatrix( 4*numElements, 1, mxREAL );
    double *F = mxGetDoubles(plhs[0]);

    for ( int el = 0; el < numElements; el++ ) {
        int vert1 = T[el]-1;
        int vert2 = T[el+numElements]-1;
        int vert3 = T[el+numElements*2]-1;
        
        int fID = 4*el;
        double DmInv1_1 = DmInv[fID];
        double DmInv2_1 = DmInv[fID+1];
        double DmInv1_2 = DmInv[fID+2];
        double DmInv2_2 = DmInv[fID+3];

        double p01 = V[vert1];
        double p02 = V[vert1 + numVerts];
        double p03 = V[vert1 + numVerts*2];
        double p11 = V[vert2];
        double p12 = V[vert2 + numVerts];
        double p13 = V[vert2 + numVerts*2];
        double p21 = V[vert3];
        double p22 = V[vert3 + numVerts];
        double p23 = V[vert3 + numVerts*2];

        double t2 = -p11;
        double t3 = -p12;
        double t4 = -p13;
        double t5 = -p21;
        double t6 = -p22;
        double t7 = -p23;
        double t8 = p01+t2;
        double t9 = p02+t3;
        double t10 = p03+t4;
        double t11 = p01+t5;
        double t12 = p02+t6;
        double t13 = p03+t7;
        double t14 = t8*t8;
        double t15 = t9*t9;
        double t16 = t10*t10;
        double t17 = t11*t11;
        double t18 = t12*t12;
        double t19 = t13*t13;
        double t20 = t8*t11;
        double t21 = t9*t12;
        double t22 = t10*t13;
        double t23 = t14+t15+t16;
        double t24 = t17+t18+t19;
        double t25 = t20+t21+t22;

        F[fID] = DmInv1_1*t23+DmInv2_1*t25;
        F[fID+1] = DmInv1_2*t23+DmInv2_2*t25;
        F[fID+2] = DmInv1_1*t25+DmInv2_1*t24;
        F[fID+3] = DmInv1_2*t25+DmInv2_2*t24;
    }

}
