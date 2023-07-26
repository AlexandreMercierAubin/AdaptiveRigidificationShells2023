function [restFrame, area] = makeElements3D(V, T)
% MAKEELEMENTS Creates elements from tetrahedrons.
%   makeElements( V, T ) returns an array of structures given 3D
%   points V (stored as rows), and tetrahedrons T defined by indices
%   (1 indexed and stored as rows).

restFrame = zeros(3,3,size(T,1));
area = zeros(size(T,1),1);
if size(T,2) == 3
    area = triangleArea3D(V,T);
end

for index = 1:size(T, 1)
    v1id = T(index, 1);
    v2id = T(index, 2);
    v3id = T(index, 3);
    if size(T,2) == 4
        v4id = T(index, 4);
        restFrame(:,:,index) = [V(v2id, :)',V(v3id, :)',V(v4id, :)'] - V(v1id, :)';
        area(index) = abs(det(restFrame(:,:,index)) / 6);
    else
        x = V(v2id, :)' - V(v1id, :)';
        y = V(v3id, :)' - V(v1id, :)';
        restFrame(:,:,index) = [x, y, zeros(3,1)];
%         n = cross(x,y);
    end
end