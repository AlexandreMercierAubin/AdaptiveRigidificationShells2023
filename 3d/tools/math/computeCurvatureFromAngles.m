function [curvature, angles, DualEdgeLength] = computeCurvatureFromAngles(meshTri3D, p,dualEdgeType, bendType,angles)
%COMPUTECURVATURE [curvature, angles] = computeCurvatureFromAngles(meshTri3D, p)
    % computes the curvature of a mesh    
    if nargin <3
        dualEdgeType = 1;
    end
    if nargin < 4
        bendType = 1;
    end
    
    V = meshTri3D.formatPositions(p);
    T = meshTri3D.t(:,1:3);
    BC = mexBarycenters(V, T);

    euclidean = 2;
    if bendType ~= 1
        curvature = angles;
        if bendType == 2
            return;
        end
    end

    if dualEdgeType == 1
        DualEdgeLength = vecnorm(BC(meshTri3D.bendingEdgesElements(:,1),:)-BC(meshTri3D.bendingEdgesElements(:,2),:), euclidean,2);
    elseif dualEdgeType == 2 % This second version is much more stable when static, but not when moving
        DualEdgeLength = vecnorm(BC(meshTri3D.bendingEdgesElements(:,1),:),euclidean,2) + vecnorm(BC(meshTri3D.bendingEdgesElements(:,2),:),euclidean,2);
    elseif dualEdgeType == 3
        DualEdgeLength = meshTri3D.restEdgeAverageHeight;
    else
        DualEdgeLength1 = vecnorm(BC(meshTri3D.bendingEdgesElements(:,1),:)-BC(meshTri3D.bendingEdgesElements(:,2),:), euclidean,2);
        DualEdgeLength2 = vecnorm(BC(meshTri3D.bendingEdgesElements(:,1),:),euclidean,2) + vecnorm(BC(meshTri3D.bendingEdgesElements(:,2),:),euclidean,2);
        DualEdgeLength = min([DualEdgeLength1,DualEdgeLength2],[],2);
    end
    if bendType == 1
        curvature = 2*angles./DualEdgeLength;
    end
end
