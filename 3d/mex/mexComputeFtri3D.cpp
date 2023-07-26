/* 
 * [ N ] = mexComputeFtri3D( V, T , DmInv);
 *
 * Input:
 *  V        #vertices by 3 matrix of vertex position
 *  T        #E by 3 vertices of elements
 *  N        deformed normals
  *  D        dphidx entries 
 * 
 * Output:
 *  F         #F 4el vector of 2 by 2 deformation gradients
 *
 * To compile type: mex -R2018a mexComputeFtri3D.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h> 

/* 
 * The gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ) {
        mexErrMsgIdAndTxt("ARP:mexComputeFtri3D:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    
    const size_t numDOFs = mxGetM(prhs[0]);
    double *V = mxGetDoubles(prhs[0]);
    const size_t numVerts = mxGetM(prhs[0]);
    int* T = (int*)mxGetData(prhs[1]);       //#tri by 3 indices (starting from 1) of verts of each element
    const size_t numElements = mxGetM(prhs[1]);  
    double* N = mxGetDoubles(prhs[2]);
    const size_t numNorms = mxGetM(prhs[2]);
    double* D = mxGetDoubles(prhs[3]);

    if ( numNorms != numElements ) {
        mexErrMsgIdAndTxt("ARP:mexComputeFtri3D:normalsElems","There should be as many normals as elements");
    }

    // create the output vectors 
    plhs[0] = mxCreateDoubleMatrix( 9*numElements, 1, mxREAL );
    double *F = mxGetDoubles(plhs[0]);

    for ( int el = 0; el < numElements; el++ ) {
        int vert1 = T[el]-1;
        int vert2 = T[el+numElements]-1;
        int vert3 = T[el+numElements*2]-1;
        
        int fID = 9*el;
		int dID = 12*el;
        double D1_1 = D[dID];
        double D2_1 = D[dID+1];
        double D3_1 = D[dID+2];
        double D4_1 = D[dID+3];
		double D1_2 = D[dID+4];
		double D2_2 = D[dID+5];
		double D3_2 = D[dID+6];
		double D4_2 = D[dID+7];
		double D1_3 = D[dID+8];
		double D2_3 = D[dID+9];
		double D3_3 = D[dID+10];
		double D4_3 = D[dID+11];

        double p01 = V[vert1];
        double p02 = V[vert1 + numVerts];
        double p03 = V[vert1 + numVerts*2];
        double p11 = V[vert2];
        double p12 = V[vert2 + numVerts];
        double p13 = V[vert2 + numVerts*2];
        double p21 = V[vert3];
        double p22 = V[vert3 + numVerts];
        double p23 = V[vert3 + numVerts*2];
		
		double n1 = N[el];
		double n2 = N[el + numElements];
		double n3 = N[el + numElements *2];

        F[fID]   = D4_1*n1+D1_1*p01+D2_1*p11+D3_1*p21;
		F[fID+1] = D4_1*n2+D1_1*p02+D2_1*p12+D3_1*p22;
		F[fID+2] = D4_1*n3+D1_1*p03+D2_1*p13+D3_1*p23;
		F[fID+3] = D4_2*n1+D1_2*p01+D2_2*p11+D3_2*p21;
		F[fID+4] = D4_2*n2+D1_2*p02+D2_2*p12+D3_2*p22;
		F[fID+5] = D4_2*n3+D1_2*p03+D2_2*p13+D3_2*p23;
		F[fID+6] = D4_3*n1+D1_3*p01+D2_3*p11+D3_3*p21;
		F[fID+7] = D4_3*n2+D1_3*p02+D2_3*p12+D3_3*p22;
		F[fID+8] = D4_3*n3+D1_3*p03+D2_3*p13+D3_3*p23;
    }

}
