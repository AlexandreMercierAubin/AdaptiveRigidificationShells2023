function [angle] = symDihedralAngle(v1,v2,o1,o2)
%DIHEDRALANGLE compute the dihedralAngle from two triangles separated by an
%edge. v1 & v2 vertices of edge. o1 & o2 vertices opposite of edge 
    b1 = o1 - v1;
    b2 = v2 - v1;
    b3 = o2 - v1;
    n1 = fastCross(b1,b2);
    n2 = fastCross(b2,b3);
    
    normn1 = fastVectorNorm(n1);
    normn2 = fastVectorNorm(n2);

    n1 = n1./normn1;
    n2 = n2./normn2;
    
    y = fastCross(n1,n2)/(normn1*normn2);

    angle = asin(y);
end

