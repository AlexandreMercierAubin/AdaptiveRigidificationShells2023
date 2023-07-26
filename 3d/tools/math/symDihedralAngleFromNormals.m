function [angle] = symDihedralAngleFromNormals(n1,n2,nCenter)
%DIHEDRALANGLE compute the dihedralAngle from two normals separated by an
%edge. 
    m1 = fastCross(n1,n2);
    y = fastDot(m1,nCenter);
    x = fastDot(n1,n2);
    angle = atan2(y,x);
%     angle = acos(n1'*n2/norm(n1)*norm(n2));
%     angle = asin(norm(fastCross(n1,n2))/norm(n1)*norm(n2));
end

