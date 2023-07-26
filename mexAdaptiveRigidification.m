% Compile the cpp mex code needed for 2D and 3D adaptive rigidificaiton
rootPwd = pwd;
cd 3d/mex
try
    mex -R2018a mexComputeSTVKGradHess2D.cpp
    mex -R2018a mexEdotNorm2D.cpp
	mex -R2018a mexEdiffNorm2D.cpp
    mex -R2018a mexFindMatches2D.cpp
    mex -R2018a mexPGS2D.cpp
    mex -R2018a mexPGS2DwithJAinvJT.cpp
    mex -R2018a mexRigidBodyProperties2D.cpp
    mex -R2018a mexRigidConnectedComponents2D.cpp
    mex -R2018a mexNeoHookean2D.cpp
    mex -R2018a mexSTVK3D.cpp
    mex -R2018a mexNeoHookean3D.cpp
    mex -R2018a mexCorotational3D.cpp
    mex -R2018a mexEdotNorm3D.cpp
    mex -R2018a mexFindMatches3D.cpp
    mex -R2018a mexPGS3D.cpp
    mex -R2018a mexPGS3DwithJAinvJT.cpp
    mex -R2018a mexRigidBodyProperties3D.cpp
    mex -R2018a mexRigidBodyConnectedComponents3D.cpp
	%mex -R2018a mexAngleBetweenVectors.cpp  -I"../../lib/gptoolbox/mex/external/libigl/external/eigen" -I"../../lib/gptoolbox/mex/external/libigl/include"
	mex -R2018a mexAddShellNormalDeformation.cpp
	mex -R2018a mexRigidBodyConnectedMixed.cpp
    mex -R2018a mexGrinspunPlanar.cpp
    mex -R2018a mexNormals.cpp
    mex -R2018a mexBarycenters.cpp
    mex -R2018a mexGrinspunBendingGradHess.cpp
    mex -R2018a mexComputeFtri3D.cpp
    %uncomment this one for efficient shells strain limiting
%     mex -R2018a mexStrainLimiting.cpp -I"../../lib/gptoolbox/mex/external/libigl/external/eigen" -I"../../lib/gptoolbox/mex/external/libigl/include"
catch
    cd ../..    
end
cd ../..
cd util/mex
try
    mex -R2018a mMD5.c
catch
    cd ../..    
end
cd ../..
%-------you can temporarily comment this if you haven't built ipc-toolkit
%ipc or ECCD won't work
% ipctoolkitIncludePath = ['-I' fullfile(rootPwd, 'lib','ipc-toolkit','build','include')];
% ipctoolkitLibPath = ['-L' fullfile(rootPwd, 'lib','ipc-toolkit','build','Release')];
% ipctoolkitlibPath = ['-l','IPCToolkit.lib'];
% eigenpath = ['-I' fullfile(rootPwd, 'lib','eigen')];
% mex('-v', '-R2018a', eigenpath, ipctoolkitIncludePath, 'pointTriangleCCD.cpp', ipctoolkitLibPath, ipctoolkitlibPath);
% mex('-v', '-R2018a', eigenpath, ipctoolkitIncludePath, 'mexIPC.cpp', ipctoolkitLibPath, ipctoolkitlibPath);

%-------

%cd ../..
