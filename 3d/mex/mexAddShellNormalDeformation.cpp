/*==========================================================
 * computes the angle between two vectors wrt a rotation axis
 *
 * [ ii, jj , dndpVal, valCount , n_vec] = mexAddShellNormalDeformation( V, C, isElementShell, materialSpaceNormals );

 * To compile type: mex -R2018a mexAddShellNormalDeformation.cpp
 *========================================================*/

#include "mex.h"
#include <math.h>

/* 
 * The gateway function
 */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    
    if ( nrhs != 4 ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:nrhs","4 inputs required.");
    }
    // make sure the arguments are dense real double
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    if ( mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustBeDenseRealDouble","Parameters must be dense int.");
    }
    if ( !mxIsLogical(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustBeLogical","Parameter isElementShell must be logical");
    }
    if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustBeDenseRealDouble","Parameters must be dense real double.");
    }
    if ( mxGetN(prhs[0]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustHave3cols","Parameter V must have 3 columns.");
    }
    if ( mxGetN(prhs[1]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustHave3cols","Parameter T must have 3 columns.");
    }
    if ( mxGetN(prhs[3]) != 3 ) {
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:mustHave3cols","Parameter materialSpaceNormals must have 3 columns.");
    }
    if (mxGetM(prhs[1]) != mxGetM(prhs[2])|| mxGetM(prhs[2]) != mxGetM(prhs[3])){
        mexErrMsgIdAndTxt("ARP:mexAddShellNormalDeformation:equalRows","Parameter T,isShellElement,and materialSpaceNormals must have the same number of rows.");
    }

    double *V = mxGetDoubles(prhs[0]);
    const size_t nvertices = mxGetM(prhs[0]);
    int* T = (int*)mxGetData(prhs[1]);
    const size_t nelements = mxGetM(prhs[1]);
    bool* isElementShell = mxGetLogicals(prhs[2]);
    double *referenceSpaceNormals = mxGetDoubles(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix( 81*nelements, 1, mxREAL );
    double *ii = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 81*nelements, 1, mxREAL );
    double *jj = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 81*nelements, 1, mxREAL );
    double *vals = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 1, 1, mxREAL );
    double *valCount = mxGetDoubles(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix( nelements*9, 1, mxREAL );
    double *n_vec = mxGetDoubles(plhs[4]);
    int c = 0;

    for(size_t i = 0; i < nelements; ++i) {
        if(!isElementShell[i]){
            n_vec[i*9] = 0; n_vec[i*9+1] = 0;n_vec[i*9+2] = 0;n_vec[i*9+3] = 0;n_vec[i*9+4] = 0;n_vec[i*9+5] = 0;n_vec[i*9+6] = 0;n_vec[i*9+7] = 0;n_vec[i*9+8] = 0;
            continue;
        }
        int v1ID = T[i];
        int v2ID = T[i+nelements];
        int v3ID = T[i+2*nelements];
        double v11 = V[(v1ID-1)];
        double v12 = V[(v1ID-1)+nvertices];
        double v13 = V[(v1ID-1)+(nvertices*2)];

        double v21 = V[(v2ID-1)];
        double v22 = V[(v2ID-1)+nvertices];
        double v23 = V[(v2ID-1)+(nvertices*2)];

        double v31 = V[(v3ID-1)];
        double v32 = V[(v3ID-1)+nvertices];
        double v33 = V[(v3ID-1)+(nvertices*2)];

        double n1 = referenceSpaceNormals[i];
        double n2 = referenceSpaceNormals[i+nelements];
        double n3 = referenceSpaceNormals[i+(2*nelements)];

        double t2 = -v21;
        double t3 = -v22;
        double t4 = -v23;
        double t5 = -v31;
        double t6 = -v32;
        double t7 = -v33;
        double t8 = t2+v11;
        double t9 = t3+v12;
        double t10 = t4+v13;
        double t11 = t5+v11;
        double t12 = t6+v12;
        double t13 = t7+v13;
        double t14 = t5+v21;
        double t15 = t6+v22;
        double t16 = t7+v23;
        double t17 = t8*t12;
        double t18 = t9*t11;
        double t19 = t8*t13;
        double t20 = t10*t11;
        double t21 = t9*t13;
        double t22 = t10*t12;
        double t23 = -t18;
        double t24 = -t20;
        double t25 = -t22;
        double t26 = t17+t23;
        double t27 = t19+t24;
        double t28 = t21+t25;
        double t29 = fabs(t26);
        double t30 = fabs(t27);
        double t31 = fabs(t28);
        double t32 = t26*t26;
        double t33 = t27*t27;
        double t34 = t28*t28;
        double t35 = t29*t29;
        double t36 = t30*t30;
        double t37 = t31*t31;
        double t38 = t35+t36+t37;
        double t39 = 1.0/t38;
        double t40 = 1.0/sqrt(t38);
        double t41 = t40*t40*t40;
        double t42 = t32*t39;
        double t43 = t33*t39;
        double t44 = t34*t39;
        double t45 = t42-1.0;
        double t46 = t43-1.0;
        double t47 = t44-1.0;
        double t48 = t8*t26*t27*t41;
        double t49 = t9*t26*t27*t41;
        double t50 = t8*t26*t28*t41;
        double t51 = t10*t26*t27*t41;
        double t52 = t9*t26*t28*t41;
        double t53 = t8*t27*t28*t41;
        double t54 = t10*t26*t28*t41;
        double t55 = t9*t27*t28*t41;
        double t56 = t11*t26*t27*t41;
        double t57 = t10*t27*t28*t41;
        double t58 = t12*t26*t27*t41;
        double t59 = t11*t26*t28*t41;
        double t60 = t13*t26*t27*t41;
        double t61 = t12*t26*t28*t41;
        double t62 = t11*t27*t28*t41;
        double t63 = t13*t26*t28*t41;
        double t64 = t12*t27*t28*t41;
        double t65 = t14*t26*t27*t41;
        double t66 = t13*t27*t28*t41;
        double t67 = t15*t26*t27*t41;
        double t68 = t14*t26*t28*t41;
        double t69 = t16*t26*t27*t41;
        double t70 = t15*t26*t28*t41;
        double t71 = t14*t27*t28*t41;
        double t72 = t16*t26*t28*t41;
        double t73 = t15*t27*t28*t41;
        double t74 = t16*t27*t28*t41;
        double t75 = -t50;
        double t76 = -t54;
        double t77 = -t57;
        double t78 = -t59;
        double t79 = -t63;
        double t80 = -t66;
        double t81 = -t68;
        double t82 = -t72;
        double t83 = -t74;
        double t84 = t8*t40*t45;
        double t85 = t9*t40*t45;
        double t86 = t8*t40*t46;
        double t87 = t10*t40*t46;
        double t88 = t11*t40*t45;
        double t89 = t9*t40*t47;
        double t90 = t12*t40*t45;
        double t91 = t10*t40*t47;
        double t92 = t11*t40*t46;
        double t93 = t13*t40*t46;
        double t94 = t14*t40*t45;
        double t95 = t12*t40*t47;
        double t96 = t15*t40*t45;
        double t97 = t13*t40*t47;
        double t98 = t14*t40*t46;
        double t99 = t16*t40*t46;
        double t100 = t15*t40*t47;
        double t101 = t16*t40*t47;
        double t102 = t48+t52;
        double t103 = t52+t57;
        double t104 = t56+t61;
        double t105 = t61+t66;
        double t106 = t65+t70;
        double t107 = t70+t74;
        double t108 = t48+t77;
        double t109 = t56+t80;
        double t110 = t65+t83;
        double t111 = t51+t85;
        double t112 = t49+t87;
        double t113 = t55+t86;
        double t114 = t53+t89;
        double t115 = t60+t90;
        double t116 = t58+t93;
        double t117 = t64+t92;
        double t118 = t62+t95;
        double t119 = t69+t96;
        double t120 = t67+t99;
        double t121 = t73+t98;
        double t122 = t71+t100;
        double t123 = t76+t84;
        double t124 = t75+t91;
        double t125 = t79+t88;
        double t126 = t78+t97;
        double t127 = t82+t94;
        double t128 = t81+t101;

        //entry 0
        vals[c] = -n1*t107;
        ii[c] = i*9+1;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = -n2*t107;
        ii[c] = i*9+1+1;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = -n3*t107;
        ii[c] = i*9+1+2;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = n1*t120;
        ii[c] = i*9+1+3;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = n2*t120;
        ii[c] = i*9+1+4;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = n3*t120;
        ii[c] = i*9+1+5;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = -n1*t119;
        ii[c] = i*9+1+6;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = -n2*t119;
        ii[c] = i*9+1+7;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = -n3*t119;
        ii[c] = i*9+1+8;
        jj[c] = v1ID*3-2;
        c++;
        vals[c] = n1*(t68-t101);
        ii[c] = i*9+1;
        jj[c] = v1ID*3-1;
        c++;
        //entry 10
        vals[c] = n2*(t68-t101);
        ii[c] = i*9+1+1;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = n3*(t68-t101);
        ii[c] = i*9+1+2;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = -n1*t110;
        ii[c] = i*9+1+3;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = -n2*t110;
        ii[c] = i*9+1+4;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = -n3*t110;
        ii[c] = i*9+1+5;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = -n1*(t72-t94);
        ii[c] = i*9+1+6;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = -n2*(t72-t94);
        ii[c] = i*9+1+7;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = -n3*(t72-t94);
        ii[c] = i*9+1+8;
        jj[c] = v1ID*3-1;
        c++;
        vals[c] = n1*t122;
        ii[c] = i*9+1;
        jj[c] = v1ID*3;
        c++;
        vals[c] = n2*t122;
        ii[c] = i*9+1+1;
        jj[c] = v1ID*3;
        c++;
        //entry 20
        vals[c] = n3*t122;
        ii[c] = i*9+1+2;
        jj[c] = v1ID*3;
        c++;
        vals[c] = -n1*t121;
        ii[c] = i*9+1+3;
        jj[c] = v1ID*3;
        c++;
        vals[c] = -n2*t121;
        ii[c] = i*9+1+4;
        jj[c] = v1ID*3;
        c++;
        vals[c] = -n3*t121;
        ii[c] = i*9+1+5;
        jj[c] = v1ID*3;
        c++;
        vals[c] = n1*t106;
        ii[c] = i*9+1+6;
        jj[c] = v1ID*3;
        c++;
        vals[c] = n2*t106;
        ii[c] = i*9+1+7;
        jj[c] = v1ID*3;
        c++;
        vals[c] = n3*t106;
        ii[c] = i*9+1+8;
        jj[c] = v1ID*3;
        c++;
        vals[c] = n1*t105;
        ii[c] = i*9+1;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = n2*t105;
        ii[c] = i*9+1+1;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = n3*t105;
        ii[c] = i*9+1+2;
        jj[c] = v2ID*3-2;
        c++;
        //entry 30
        vals[c] = -n1*t116;
        ii[c] = i*9+1+3;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = -n2*t116;
        ii[c] = i*9+1+4;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = -n3*t116;
        ii[c] = i*9+1+5;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = n1*t115;
        ii[c] = i*9+1+6;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = n2*t115;
        ii[c] = i*9+1+7;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = n3*t115;
        ii[c] = i*9+1+8;
        jj[c] = v2ID*3-2;
        c++;
        vals[c] = -n1*(t59-t97);
        ii[c] = i*9+1;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = -n2*(t59-t97);
        ii[c] = i*9+1+1;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = -n3*(t59-t97);
        ii[c] = i*9+1+2;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = n1*t109;
        ii[c] = i*9+1+3;
        jj[c] = v2ID*3-1;
        c++;
        //entry 40
        vals[c] = n2*t109;
        ii[c] = i*9+1+4;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = n3*t109;
        ii[c] = i*9+1+5;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = n1*(t63-t88);
        ii[c] = i*9+1+6;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = n2*(t63-t88);
        ii[c] = i*9+1+7;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = n3*(t63-t88);
        ii[c] = i*9+1+8;
        jj[c] = v2ID*3-1;
        c++;
        vals[c] = -n1*t118;
        ii[c] = i*9+1;
        jj[c] = v2ID*3;
        c++;
        vals[c] = -n2*t118;
        ii[c] = i*9+1+1;
        jj[c] = v2ID*3;
        c++;
        vals[c] = -n3*t118;
        ii[c] = i*9+1+2;
        jj[c] = v2ID*3;
        c++;
        vals[c] = n1*t117;
        ii[c] = i*9+1+3;
        jj[c] = v2ID*3;
        c++;
        vals[c] = n2*t117;
        ii[c] = i*9+1+4;
        jj[c] = v2ID*3;
        c++;
        //entry 50
        vals[c] = n3*t117;
        ii[c] = i*9+1+5;
        jj[c] = v2ID*3;
        c++;
        vals[c] = -n1*t104;
        ii[c] = i*9+1+6;
        jj[c] = v2ID*3;
        c++;
        vals[c] = -n2*t104;
        ii[c] = i*9+1+7;
        jj[c] = v2ID*3;
        c++;
        vals[c] = -n3*t104;
        ii[c] = i*9+1+8;
        jj[c] = v2ID*3;
        c++;
        vals[c] = -n1*t103;
        ii[c] = i*9+1;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = -n2*t103;
        ii[c] = i*9+1+1;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = -n3*t103;
        ii[c] = i*9+1+2;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = n1*t112;
        ii[c] = i*9+1+3;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = n2*t112;
        ii[c] = i*9+1+4;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = n3*t112;
        ii[c] = i*9+1+5;
        jj[c] = v3ID*3-2;
        c++;
        //entry 60
        vals[c] = -n1*t111;
        ii[c] = i*9+1+6;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = -n2*t111;
        ii[c] = i*9+1+7;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = -n3*t111;
        ii[c] = i*9+1+8;
        jj[c] = v3ID*3-2;
        c++;
        vals[c] = n1*(t50-t91);
        ii[c] = i*9+1;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = n2*(t50-t91);
        ii[c] = i*9+1+1;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = n3*(t50-t91);
        ii[c] = i*9+1+2;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = -n1*t108;
        ii[c] = i*9+1+3;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = -n2*t108;
        ii[c] = i*9+1+4;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = -n3*t108;
        ii[c] = i*9+1+5;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = -n1*(t54-t84);
        ii[c] = i*9+1+6;
        jj[c] = v3ID*3-1;
        c++;
        //entry 70
        vals[c] = -n2*(t54-t84);
        ii[c] = i*9+1+7;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = -n3*(t54-t84);
        ii[c] = i*9+1+8;
        jj[c] = v3ID*3-1;
        c++;
        vals[c] = n1*t114;
        ii[c] = i*9+1;
        jj[c] = v3ID*3;
        c++;
        vals[c] = n2*t114;
        ii[c] = i*9+1+1;
        jj[c] = v3ID*3;
        c++;
        vals[c] = n3*t114;
        ii[c] = i*9+1+2;
        jj[c] = v3ID*3;
        c++;
        vals[c] = -n1*t113;
        ii[c] = i*9+1+3;
        jj[c] = v3ID*3;
        c++;
        vals[c] = -n2*t113;
        ii[c] = i*9+1+4;
        jj[c] = v3ID*3;
        c++;
        vals[c] = -n3*t113;
        ii[c] = i*9+1+5;
        jj[c] = v3ID*3;
        c++;
        vals[c] = n1*t102;
        ii[c] = i*9+1+6;
        jj[c] = v3ID*3;
        c++;
        vals[c] = n2*t102;
        ii[c] = i*9+1+7;
        jj[c] = v3ID*3;
        c++;
        //entry 80
        vals[c] = n3*t102;
        ii[c] = i*9+1+8;
        jj[c] = v3ID*3;
        c++;

        n_vec[i*9] = n1*t28*t40;
        n_vec[i*9+1] = n2*t28*t40;
        n_vec[i*9+2] = n3*t28*t40;
        n_vec[i*9+3] = -n1*t27*t40;
        n_vec[i*9+4] = -n2*t27*t40;
        n_vec[i*9+5] = -n3*t27*t40;
        n_vec[i*9+6] = n1*t26*t40;
        n_vec[i*9+7] = n2*t26*t40;
        n_vec[i*9+8] = n3*t26*t40;
    }	
    valCount[0] = 1.0*c;
}
