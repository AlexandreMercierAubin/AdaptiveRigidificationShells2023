function m = computeMass(mesh)
    % computeMass computes the diagonal of the mass matrix
    %   computeMass( mesh, rho ) returns the lumped node diagonal of the mass 
    %   matrix for the given mesh

    m = zeros(mesh.N * 3, 1);
    for i = 1:size(mesh.t, 1)
        rho = mesh.materials(mesh.materialIndex(i)).rho;
        if mesh.elementType(i) == mesh.elementTypeEnum.Shell
            T = mesh.t(i,1:3);
        else
            T = mesh.t(i,:);
        end
        dofInds =[T*3-2;T*3-1;T*3];
        nverts = numel(T);
        % number of verts in the element (tri3D or tet)
        %mass is the density multiplied by the area divided by the
        %number of nodes of the element
        m(dofInds) = m(dofInds) +  (mesh.elV(i) * rho)/nverts;
    end
end