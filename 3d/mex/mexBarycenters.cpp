/* 
 * [ N ] = mexBarycenters( V, T );
 *
 * Input:
 *  V        #vertices by 3 matrix of vertex position
 *  T        #E by 3 vertices of elements
 * 
 * Output:
 *  B         #E by 3 stack of triangle barycenters
 *
 * To compile type: mex -R2018a mexBarycenters.cpp
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h> 

void addEntry(double *vals,double *valsD, double alpha1, int &counter, double * ii,double * jj, int key1, int key2, double value){
    ii[counter] = key1;
    jj[counter] = key2;
    vals[counter] = value;
    valsD[counter] = alpha1*value;
    ++counter;
}

/* 
 * The gateway function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    const size_t numDOFs = mxGetM(prhs[0]);
    double *V = mxGetDoubles(prhs[0]);
    const size_t numVerts = mxGetM(prhs[0]);
    int* T = (int*)mxGetData(prhs[1]);       //#tri by 3 indices (starting from 1) of verts of each element
    const size_t numElements = mxGetM(prhs[1]);

    // create the output vectors 
    plhs[0] = mxCreateDoubleMatrix( numElements, 3, mxREAL );
    double *B = mxGetDoubles(plhs[0]);

    for ( int el = 0; el < numElements; el++ ) {
        int vert1 = T[el]-1;
        int vert2 = T[el+numElements]-1;
        int vert3 = T[el+numElements*2]-1;
        double V1_1 = V[vert1];
        double V1_2 = V[vert1 + numVerts];
        double V1_3 = V[vert1 + numVerts*2];
        double V2_1 = V[vert2];
        double V2_2 = V[vert2 + numVerts];
        double V2_3 = V[vert2 + numVerts*2];
        double V3_1 = V[vert3];
        double V3_2 = V[vert3 + numVerts];
        double V3_3 = V[vert3 + numVerts*2];

        B[el] = V1_1+V2_1+V3_1;
        B[el+numElements] = V1_2+V2_2+V3_2;
        B[el+numElements*2] = V1_3+V2_3+V3_3;
    }

}
