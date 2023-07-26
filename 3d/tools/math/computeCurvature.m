function [curvature, angles, BC, DualEdgeLength, bendingEdgeLength] = computeCurvature(meshTri3D, p,dualEdgeType, bendType)
    if nargin <3
        dualEdgeType = 1;
    end
    if nargin < 4
        bendType = 1;
    end
    %COMPUTECURVATURE [curvature, angles, BC, DualEdgeLength, bendingEdgeLength] = computeCurvature(meshTri3D, p)
    % computes the curvature of a mesh
    V = meshTri3D.formatPositions(p);
    T = meshTri3D.t(:,1:3);
    BC = mexBarycenters(V, T);
%     BC = barycenters(V, T);
    
%      triNormals = normals(V, T, 'Stable', true);
%      triNormals = normr(triNormals);
    triNormals = mexNormals(V,T);
    
    norm1 = triNormals(meshTri3D.bendingEdgesElements(:,1),:);
    norm2 = triNormals(meshTri3D.bendingEdgesElements(:,2),:);
    vertex1full = meshTri3D.edgeNormalConsistentVertices(:,1);
    vertex2full = meshTri3D.edgeNormalConsistentVertices(:,2);
    v1full = p([vertex1full*3-2,vertex1full*3-1,vertex1full*3]);
    v2full = p([vertex2full*3-2,vertex2full*3-1,vertex2full*3]);
    % needs the reshape cause sometimes it gives column vectors when
    % there is just one edge
    v1full = reshape(v1full,[],3);
    v2full = reshape(v2full,[],3);
    
    euclidean = 2;
    e0full = v2full-v1full;
    bendingEdgeLength = vecnorm(e0full,euclidean,2);
    edgeRestNormalized = e0full./bendingEdgeLength;
    angles = mexAngleBetweenVectors(norm1,norm2,edgeRestNormalized);
%     [anglesTest] = dihedralAngleFromNormals(norm1,norm2,edgeRestNormalized);
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
