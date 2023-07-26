function [I] = computeMassRhatSquared(x,y,z,m)
    mxy = m*x*y;
    mxz = m*x*z;
    myz = m*y*z;
    mxx = m*x*-x;
    myy = m*y*-y;
    mzz = m*z*-z;
    I = [mzz+myy,   mxy,              mxz;
          mxy,                 mzz+mxx,    myz;
          mxz,                 myz,              myy+mxx];
end