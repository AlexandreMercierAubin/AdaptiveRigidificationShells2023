function computeActiveDOFs(mesh)
% active dofs will also have rigid dofs removed (override by adaptive mesh)
    mesh.activeDOFs = mesh.unpinnedDOFs;
    mesh.ActiveDofsCorrespondingID = false(mesh.N*3,1);
    mesh.ActiveDofsCorrespondingID(mesh.activeDOFs) = true;
end