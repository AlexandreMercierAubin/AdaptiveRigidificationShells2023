function [angles] = dihedralAngleFromNormals(N1,N2,EV)
%DIHEDRALANGLE compute the dihedralAngle from two normals separated by an
%edge. 
    numAngles = size(N1,1);
    angles = zeros(numAngles,1);
    for i = 1:numAngles
        NCross = fastCross(EV(i,:),N1(i,:));
        NCrossNorm = NCross./fastVectorNorm(NCross);
        E = [EV(i,:)', N1(i,:)', NCrossNorm];
        c = E'*N2(i,:)';
        angles(i) = atan2(c(3),c(2));
    end
end

