function [newB, newF]=addShellNormalDeformation(meshes, p, normals, settings)
    if settings.addShellNormalDeformation
        isShellElement = meshes.elementType == meshes.elementTypeEnum.Shell;

        V = meshes.formatPositions(p);
        %TODO: regen this mex
%         [ii,jj,vals, valCount, n_vec] = mexAddShellNormalDeformation(V, meshes.t(:,1:3), isShellElement, meshes.referenceSpaceNormals);
        [n_vec,ii,jj,vals, valCount] = addShellNormalDeformationHelper(V, meshes.t(:,1:3), isShellElement, meshes.referenceSpaceNormals, normals);

        %scales back the vector to the actual number of entries
        ii = ii(1:valCount);
        jj = jj(1:valCount);
        
        vals = vals(1:valCount);

        fullNdndp = sparse(ii,jj,vals,9*size(meshes.t,1),numel(p));
        newB = meshes.BnoNormal + fullNdndp;
        newF = meshes.BnoNormal*p + n_vec; %or newB here?
    else
        isShellElement = meshes.elementType == meshes.elementTypeEnum.Shell;
        V = meshes.getPositionFormatted;
        T = meshes.t(isShellElement,:);
        N = normals(isShellElement,:);
        dphidx = meshes.dphidx(:,:,isShellElement);
        Ftri = mexComputeFtri3D(V,T,N,dphidx(:));
        newF = meshes.B*p;
        newF(meshes.TriFRows) = Ftri;
        newB = meshes.B;
    end
end