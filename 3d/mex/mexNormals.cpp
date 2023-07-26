/* 
 * [ N ] = mexNormals( V, T );
 *
 * Input:
 *  V        #vertices by 3 matrix of vertex position
 *  T        #E by 3 vertices of elements
 * 
 * Output:
 *  N         #E by 3 stack of normals
 *
 * To compile type: mex -R2018a mexNormals.cpp
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
    double *N = mxGetDoubles(plhs[0]);

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

        // code gen
                double t2 = -V2_1;
        double t3 = -V2_2;
        double t4 = -V2_3;
        double t5 = -V3_1;
        double t6 = -V3_2;
        double t7 = -V3_3;
        double t8 = V1_1+t2;
        double t9 = V1_2+t3;
        double t10 = V1_3+t4;
        double t11 = V1_1+t5;
        double t12 = V1_2+t6;
        double t13 = V1_3+t7;
        double t14 = V2_1+t5;
        double t15 = V2_2+t6;
        double t16 = V2_3+t7;
        double t17 = (t8*t12)/3.0;
        double t18 = (t9*t11)/3.0;
        double t19 = (t8*t13)/3.0;
        double t20 = (t10*t11)/3.0;
        double t21 = (t9*t13)/3.0;
        double t22 = (t10*t12)/3.0;
        double t23 = (t8*t15)/3.0;
        double t24 = (t9*t14)/3.0;
        double t25 = (t8*t16)/3.0;
        double t26 = (t10*t14)/3.0;
        double t27 = (t9*t16)/3.0;
        double t28 = (t10*t15)/3.0;
        double t29 = (t11*t15)/3.0;
        double t30 = (t12*t14)/3.0;
        double t31 = (t11*t16)/3.0;
        double t32 = (t13*t14)/3.0;
        double t33 = (t12*t16)/3.0;
        double t34 = (t13*t15)/3.0;
        double t35 = -t18;
        double t36 = -t20;
        double t37 = -t22;
        double t38 = -t24;
        double t39 = -t26;
        double t40 = -t28;
        double t41 = -t30;
        double t42 = -t32;
        double t43 = -t34;
        double t44 = t17+t23+t29+t35+t38+t41;
        double t45 = t19+t25+t31+t36+t39+t42;
        double t46 = t21+t27+t33+t37+t40+t43;
        double t47 = fabs(t44);
        double t48 = fabs(t45);
        double t49 = fabs(t46);
        double t50 = t47*t47;
        double t51 = t48*t48;
        double t52 = t49*t49;
        double t53 = t50+t51+t52;
        double t54 = 1.0/sqrt(t53);
        N[el]  = t46*t54;
        N[el+numElements] = -t45*t54;
        N[el+numElements*2] = t44*t54;
    }

}
