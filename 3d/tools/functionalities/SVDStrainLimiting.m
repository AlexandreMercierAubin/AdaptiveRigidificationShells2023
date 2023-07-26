function [newp]=SVDStrainLimiting(meshes, p, cache, settings)
    %SVDStrainLimiting(meshes, cache, settings)
    %using material properties, limits strain by SVD decomposition and
    %analyzing eigenvalues
    
    %TODO: this can cause pinned dofs to move due to strain limits. Make
    %sure that pinned dofs don't move
    newp = p;
    if settings.StrainLimitingEnabled == 2
        if isa(meshes,"AdaptiveMesh3D") %skip rigid elements as they are non-deforming
            elasticElementIDs = find(meshes.isTetElastic);
        else
            elasticElementIDs = 1:size(meshes.t,1);
        end
        [newF,newp] = mexStrainLimiting(cache.F,meshes.elUpper,meshes.elLower,meshes.mass,meshes.t,meshes.p,meshes.restFrame(:),elasticElementIDs);
        cache.F = newF;
    elseif settings.StrainLimitingEnabled == 1
        %based on http://graphics.berkeley.edu/papers/Wang-MRI-2010-12/Wang-MRI-2010-12.pdf
        numt = size(meshes.t,1);
        reshapedF = reshape(cache.F,3,3,[]);
        if isa(meshes,"AdaptiveMesh3D") %skip rigid elements as they are non-deforming
            elasticElementIDs = find(meshes.isTetElastic);
        else
            elasticElementIDs = 1:size(meshes.t,1);
        end

        for idIndex = 1:numel(elasticElementIDs)
            i = elasticElementIDs(idIndex);
            localF = reshapedF(:,:,i);
            [U,Sdiag,V] = svd(localF,"vector");
            upperbound = meshes.elUpper(i);
            lowerbound = meshes.elLower(i);

            %clamp the eigenvalues of F
            SStarDiag = Sdiag;
            SStarDiag(SStarDiag > upperbound) = upperbound;
            SStarDiag(SStarDiag < lowerbound) = lowerbound;

            %strain limiting inversion safeguard
            %should this be before clamping?
            if det(localF) < 0 
                minSdiag = min(SStarDiag);
                SStarDiag(SStarDiag == minSdiag) = -minSdiag;
            end

            %reconstructing a corrected deformation gradient
            localFStar = U*diag(SStarDiag)*V';
            cache.F(i*9-8:i*9) = localFStar(:);

            triangleVert = meshes.t(i,:);
            Dr = meshes.restFrame(:,:,i);
            if meshes.elementType(i) == meshes.elementTypeEnum.Shell
                Dr = Dr(:,1:2);
                triangleVert = triangleVert(:,1:3);
            end
            pinnedVerts = meshes.pinned(triangleVert);
            pinnedDofs = logical(reshape(repmat(pinnedVerts,1,3)',[],1));
            ids = [triangleVert*3-2; triangleVert*3-1;triangleVert*3];
            currentPos = p(ids);

            Dx = localFStar*Dr;

            %center of mass
            mass = meshes.mass(triangleVert*3);%mass is repeated per dofs so 3x N particles
            massEnd = mass(2:end);

            sumMass = sum(mass);
            DxMass = Dx*massEnd;
            weightedDistance = DxMass./sumMass;

            
            weightedPos = currentPos * mass;
            centerOfMass = weightedPos ./ sumMass;
            x0 = centerOfMass - weightedDistance;

            xn = [x0,x0+Dx];
            xn(pinnedDofs) = currentPos(pinnedDofs);
            newp(ids) = xn;
        end
    elseif settings.StrainLimitingEnabled == 3 %edge strain limit
        edges = meshes.edges;
        l0 = meshes.edgesRest;
        V = meshes.formatPositions(p);
        edgeVector = V(meshes.edges(:,1),:)-V(meshes.edges(:,2),:);
        edgeCenter = (V(meshes.edges(:,1),:)+V(meshes.edges(:,2),:))./2;
        edgeLength = vecnorm(edgeVector,2,2);
        edgeNormal = edgeVector./edgeLength;
        smax = 1.001*l0;
        smin = 0.999*l0;
        tooSmall = find(edgeLength < smin);
        tooBig = find(edgeLength > smax);

        vert1ids = edges(tooSmall,1);
        vert2ids = edges(tooSmall,2);
%         if any(tooSmall) || any(tooBig)
%             tmp = 0
%         end
        newp([vert1ids*3-2,vert1ids*3-1,vert1ids*3]) = edgeCenter(tooSmall,:) + edgeNormal(tooSmall,:).*(smin(tooSmall)./2);
        newp([vert2ids*3-2,vert2ids*3-1,vert2ids*3]) = edgeCenter(tooSmall,:) - edgeNormal(tooSmall,:).*(smin(tooSmall)./2);

        vert1ids = edges(tooBig,1);
        vert2ids = edges(tooBig,2);
        newp([vert1ids*3-2,vert1ids*3-1,vert1ids*3]) = edgeCenter(tooBig,:) + edgeNormal(tooBig,:).*(smax(tooBig,:)./2);
        newp([vert2ids*3-2,vert2ids*3-1,vert2ids*3]) = edgeCenter(tooBig,:) - edgeNormal(tooBig,:).*(smax(tooBig,:)./2);
        newp(meshes.pinnedDOFs) = p(meshes.pinnedDOFs);
    end
    
end