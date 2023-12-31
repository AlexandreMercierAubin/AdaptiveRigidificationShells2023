function updateParticles(mesh, h, deltav)
    mesh.v(mesh.ActiveDofsCorrespondingID) = mesh.v(mesh.ActiveDofsCorrespondingID) + deltav;
    mesh.v(mesh.pinnedDOFs) = mesh.v(mesh.pinnedDOFs) * 0;
    mesh.prevp = mesh.p;
    mesh.p = mesh.p + h * mesh.v;
end