function precomputeAInverseDiagonalBlocks3D( cache, mesh, h, energyModel,settings)
    % precomputeAInverseDiagonalBlocks Builds the rest pose A matrix and 
    % precomputes just the diagonal blocks of the inverse matrix. 
    
    [bigB, F]=addShellNormalDeformation(mesh, mesh.p0,mesh.referenceSpaceNormals, settings);
    
    [psi, bendingForces,shellBendingH, shellBendingDampingH, dihedralAngles] = computeBendingGradHess(mesh, mesh.p0, mesh.referenceSpaceNormals, cache, settings);
            
    if settings.useGrinspunPlanarEnergy
        isShell = mesh.elementType == mesh.elementTypeEnum.Shell;
        [ ii, jj, vals, dpsidX , valsD, Wa, Wl] = mexGrinspunPlanar(mesh.p,mesh.t,mesh.elkl,mesh.elka,mesh.area,isShell, mesh.elAlpha1, mesh.grinspunEnergyEdgeRestLength, sum(isShell));
        grinspunPlanarH = sparse(ii,jj,-vals,mesh.N*3,mesh.N*3);
        grinspunPlanarHd = sparse(ii,jj,-valsD,mesh.N*3,mesh.N*3);
        grinspunPlanarForces = -dpsidX;
        allpsi = Wa + Wl;
        sumWa = sum(Wa);
        sumWl = sum(Wl);
        psi = sumWa + sumWl;
        K = grinspunPlanarH + shellBendingH;
        Kd = grinspunPlanarHd + shellBendingDampingH;
    else
        energyModel.computeEnergy(mesh, F);
        C = energyModel.derivative2HessianC;
        dpsidF = energyModel.derivative1Gradient;
        elasticForces = energyModel.elasticForces; 
        K = energyModel.stiffness + shellBendingH;
        Kd = energyModel.stiffnessDamping + shellBendingDampingH;
    end 
    
    
    M = mesh.M;
    Md = mesh.Md;
    
    A = M - h * (-Md + Kd) - h^2 * K; 

    if isempty(cache.AInvBlocks)
        % this is ONLY needed for boundary nodes... could make faster...
        % Also, at the point that the precomp is called, this cache
        % should be empty
        cache.AInvBlocks = computeAInverseDiagonalBlocks3D( A );
    else
        disp('already have a cached AinvBlocks.. why is precmop getting called?');
    end
    
end