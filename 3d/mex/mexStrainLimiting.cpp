/*==========================================================
 *
 * To compile type: mex -R2018a mexStrainLimiting.cpp -I"../../lib/gptoolbox/mex/external/libigl/external/eigen" -I"../../lib/gptoolbox/mex/external/libigl/include"
 *========================================================*/

#include "mex.h"
#include "blas.h"
#include <math.h>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *F = mxGetDoubles(prhs[0]);
    size_t sizeF = mxGetM(prhs[0]);
    double *u = mxGetDoubles(prhs[1]);
    size_t numElements = mxGetM(prhs[1]);
    double *l = mxGetDoubles(prhs[2]);
    double *mass = mxGetDoubles(prhs[3]);
    int* T = (int*)mxGetData(prhs[4]);
    double *x = mxGetDoubles(prhs[5]);
    size_t numDOFs = mxGetM(prhs[5]);
    double* Dr = mxGetDoubles(prhs[6]);
    double* elasticElementsList = mxGetDoubles(prhs[7]);
    size_t numElastic = mxGetM(prhs[7]);

    if ( nrhs != 8 ) {
        mexErrMsgIdAndTxt("ARP:mexStrainLimiting:nrhs","8 inputs required.");
    }
    // make sure the arguments are dense real double
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ) {
        mexErrMsgIdAndTxt("ARP:mexStrainLimiting:mustBeDenseRealDouble","Parameters must be dense real double.");
    }

    plhs[0] = mxCreateDoubleMatrix(sizeF,1,mxREAL);
    double *newF = mxGetDoubles(plhs[0]);
    memcpy(newF,F,sizeF*sizeof(double));
    plhs[1] = mxCreateDoubleMatrix( numDOFs, 1, mxREAL );
    double *newx = mxGetDoubles(plhs[1]);
    memcpy(newx,x,numDOFs*sizeof(double));

    Matrix3d localF;//doesn't matter what is in it as it gets overwritten
    MatrixXd localDr(3,2);
    Vector2d massEnd;
    Vector3d weightedDistance;
    Matrix3d currentPos;
    Vector3d localMass;
    Vector3d centerOfMass;
    MatrixXd localDx(3, 2);
    for(size_t i = 0; i < numElastic; ++i) {
        size_t el = elasticElementsList[i]-1;
        int pos = el*9;
        int vert1 = T[el]-1;
        int vert2 = T[el+numElements]-1;
        int vert3 = T[el+numElements*2]-1;

        int dID = 9*el;
        //assumes this is cloth, otherwise I would need the last 3 entries
        localDr.col(0) << Dr[dID],Dr[dID+1],Dr[dID+2];
        localDr.col(1) << Dr[dID+3],Dr[dID+4],Dr[dID+5];

        int dof1 = vert1*3;
        int dof2 = vert1*3+1;
        int dof3 = vert1*3+2;
        int dof4 = vert2*3;
        int dof5 = vert2*3+1;
        int dof6 = vert2*3+2;
        int dof7 = vert3*3;
        int dof8 = vert3*3+1;
        int dof9 = vert3*3+2;

        currentPos << x[dof1],x[dof4],x[dof7],
                      x[dof2],x[dof5],x[dof8],
                      x[dof3],x[dof6],x[dof9];

        localF << F[pos], F[pos+3], F[pos+6],
                  F[pos+1], F[pos+4],F[pos+7],
                  F[pos+2], F[pos+5], F[pos+8];

        JacobiSVD<Matrix3d> svd(localF, Eigen::ComputeFullU | Eigen::ComputeFullV);
        
        auto singularValues = svd.singularValues().cwiseMin(u[el]).cwiseMax(l[el]).asDiagonal();
        Matrix3d localFStar = svd.matrixU()* singularValues * svd.matrixV().adjoint();
        //assumes no pinned vertices
        localDx = localFStar*localDr;
        localMass << mass[dof2], mass[dof5], mass[dof8];
        double sumMass = localMass.sum();
        massEnd << mass[dof5], mass[dof8];
        Vector3d DxMass = localDx * massEnd;
        weightedDistance = DxMass / sumMass;
        
        Vector3d weightedPos = currentPos*localMass;
        centerOfMass = weightedPos / sumMass;
        Vector3d x0 = centerOfMass - weightedDistance;
        
        newx[dof1] = x0(0);
        newx[dof2] = x0(1);
        newx[dof3] = x0(2);
        newx[dof4] = x0(0)+ localDx.coeff(0,0);
        newx[dof5] = x0(1)+ localDx.coeff(1,0);
        newx[dof6] = x0(2)+ localDx.coeff(2,0);
        newx[dof7] = x0(0)+ localDx.coeff(0,1);
        newx[dof8] = x0(1)+ localDx.coeff(1,1);
        newx[dof9] = x0(2)+ localDx.coeff(2,1);

        newF[pos]   = localFStar.coeff(0,0);
        newF[pos+1] = localFStar.coeff(1,0);
        newF[pos+2] = localFStar.coeff(2,0);
        newF[pos+3] = localFStar.coeff(0,1);
        newF[pos+4] = localFStar.coeff(1,1);
        newF[pos+5] = localFStar.coeff(2,1);
        newF[pos+6] = localFStar.coeff(0,2);
        newF[pos+7] = localFStar.coeff(1,2);
        newF[pos+8] = localFStar.coeff(2,2);
    }
}
