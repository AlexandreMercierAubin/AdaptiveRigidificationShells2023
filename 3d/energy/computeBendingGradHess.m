function [psi, fullg, shellBendingH, dampingH, angles]=computeBendingGradHess(mesh3D, p, normals, cache, settings)
    psi = 0;
    if nargout > 1
        fullg = zeros(mesh3D.N*3,1);
	    shellBendingH = sparse(mesh3D.N*3,mesh3D.N*3);
        dampingH = sparse(mesh3D.N*3,mesh3D.N*3);
        angles = [];
    end

	if settings.addBendingEnergy
        V = mesh3D.formatPositions(p);
        innerAngles = internalangles(V,mesh3D.t); %a23 a13 a12 -> v1 v2 v3

        areas = triangleArea3D(V,mesh3D.t);

        norm1 = normals(mesh3D.bendingEdgesElements(:,1),:);
        norm2 = normals(mesh3D.bendingEdgesElements(:,2),:);
        %get the axis of rotation
        vertex1full = mesh3D.edgeNormalConsistentVertices(:,1);
        vertex2full = mesh3D.edgeNormalConsistentVertices(:,2);
        v1full = p([vertex1full*3-2,vertex1full*3-1,vertex1full*3]);
        v2full = p([vertex2full*3-2,vertex2full*3-1,vertex2full*3]);
        % needs the reshape cause sometimes it gives column vectors when
        % there is just one edge
        v1full = reshape(v1full,[],3);
        v2full = reshape(v2full,[],3);

        e0full = v2full-v1full;
        e0norm = normr(e0full);
        angles = mexAngleBetweenVectors(norm1,norm2,e0norm);
        anglesDiff = angles-mesh3D.edgeRestDihedral;
        kd = mesh3D.kdShellScalingFactor;
        psi_i = kd.*anglesDiff.*anglesDiff;
        psi = sum(psi_i);

        if nargout > 1
            numEdges =size(mesh3D.bendingEdges,1);

            numAngles = size(innerAngles,1);
            innerAngle1 = innerAngles(mesh3D.bendingEdgeInnerAngleOrder(:,1) + numAngles*(mesh3D.bendingEdgeInnerAngleOrder(:,2)-1));
            innerAngleTilde1 = innerAngles(mesh3D.bendingEdgeInnerAngleOrder(:,3) + numAngles*(mesh3D.bendingEdgeInnerAngleOrder(:,4)-1));
            innerAngle2 = innerAngles(mesh3D.bendingEdgeInnerAngleOrder(:,5) + numAngles*(mesh3D.bendingEdgeInnerAngleOrder(:,6)-1));
            innerAngleTilde2 = innerAngles(mesh3D.bendingEdgeInnerAngleOrder(:,7) + numAngles*(mesh3D.bendingEdgeInnerAngleOrder(:,8)-1));
            
            area = areas(mesh3D.bendingEdgesElements(:,1));
            areaTilde = areas(mesh3D.bendingEdgesElements(:,2));
            
            [ii,jj, vals, valsD, fullg] = mexGrinspunBendingGradHess(p, mesh3D.edgeNormalConsistentVertices, mesh3D.edgeOppositeVertices, kd, anglesDiff, norm1, norm2, area, areaTilde, innerAngle1, innerAngleTilde1, innerAngle2, innerAngleTilde2,mesh3D.bendingEdgeAlpha1);
%             [ii,jj, vals, valsD, fullg] = grinspunBendingGradHessHelper(p, mesh3D.edgeNormalConsistentVertices, mesh3D.edgeOppositeVertices, kd, anglesDiff, norm1, norm2, area, areaTilde, innerAngle1, innerAngleTilde1, innerAngle2, innerAngleTilde2,mesh3D.bendingEdgeAlpha1);

            shellBendingH = sparse(ii,jj,vals,mesh3D.N*3,mesh3D.N*3); % This works because sparse build sums redundant entries
            dampingH = sparse(ii,jj,valsD,mesh3D.N*3,mesh3D.N*3);
        end
        %-------------------------------------------------------------
        %mexed version

%         [ii,jj,Hvals, forces, angles, psi] = mexGrinspunBendingGradHess(p, mesh3D.bendingEdges,mesh3D.edgeOppositeVertices, mesh3D.kdRestEdgeAverageOppositeHeight, mesh3D.edgeRestDihedral);
%         cache.shellBendingH = sparse(ii,jj,Hvals, mesh3D.N*3, mesh3D.N*3);
%         cache.elasticForces = cache.elasticForces + forces;
    end
end

