function quickSolve3D( cache, integrator, mesh3D, h, Jc, phi, settings, animationScripter, frame, cInfo)
    % QUICKSOLVE performs the single iteration of preconditionned conjugate
    % gradient required for derigidification, does not support constraints
    if nargin < 5 || isempty(Jc)
         Jc = zeros( 0, mesh3D.N*3); % no constraints
    end
    
    if ~mesh3D.EnableAutomaticRigidification
        cache.ApproximatedDeltaV = zeros(size(mesh3D.v));   
    end

    cgTic = tic;
     
    % NOTE: OPTIMIZATIONS...
    % quicksolve computes F and STVK forces and C for *all* elements and
    % this can be reused in the solve later (though only small parts are
    % needed).  All of this is cheap and linear in #TRI but gets expensive.
    
    
    % update forces
    mesh3D.f = zeros( mesh3D.N*3, 1 );
    % Gravity (z is vertical in matlab)
    %lets use forward euler on g to separate the gravity forces
    ii = mesh3D.unpinnedDOFs;

    % We'll do the warm start here, and reuse the result later! This 
    % will let us know what is new, and what forces to put on the rhs

    [warmStartLambdas, isNewContact] =  cache.findWarmStartLambdaArray;
    contactVertices = [cInfo.vertexID];
    
    if integrator.separateQuicksolveGravity && numel(cache.gravityDv) == 0
        %TODO: This needs to be updated on mass change or when dofs are pinned
        %TODO: separating gravity from the pcg solve does not work with
        %pinned dofs; fix it
        if cache.recomputeGravityDv == 1
            cache.gravityDv = zeros(size(mesh3D.f(ii)));
            gravity = zeros( mesh3D.N*3, 1 );
            gravity(3:3:end) = integrator.Gravity;
            
            cache.gravityDv = h*gravity(ii);
            cache.recomputeGravityDv = 0;
        end
        %adds gravity in the rhs for contacting vertices
%         oldContactVertices = contactVertices(~isNewContact);
        
%         gravityDv(oldContactVertices*3)=0;
%         mesh3D.f(oldContactVertices*3) = mesh3D.mass(oldContactVertices*3) * integrator.Gravity;
    elseif ~integrator.separateQuicksolveGravity
        cache.gravityDv = zeros(size(mesh3D.f(ii)));
        mesh3D.f(3:3:end) = mesh3D.mass(3:3:end) * integrator.Gravity;
    end
    gravityDv = cache.gravityDv;
    mesh3D.f(mesh3D.pinnedDOFs) = 0;

    if settings.useGrinspunPlanarEnergy %TODO: make this work with tets
        mesh3D.f = mesh3D.f + cache.grinspunPlanarForces + integrator.CustomForce;
    else
        mesh3D.f = mesh3D.f + cache.elasticForces + integrator.CustomForce;
    end
    if (settings.addBendingEnergy)
        mesh3D.f = mesh3D.f + cache.bendingForces;
    end
    
    %add animationScript forces when scripted
    for i=1:numel(animationScripter)
        mesh3D.f = animationScripter{i}.scriptForceAnimation(mesh3D.f, frame, h);
    end
    
    rhs = h * mesh3D.f(ii);
    %params = zeros(size(mesh.t, 1), 2);

    M = mesh3D.Mii;
    Md = mesh3D.Mdii;
    v = mesh3D.v(ii) + gravityDv;
    
    K = cache.K(ii,ii);
    Kd = cache.Kd(ii,ii);
    D = cache.D(ii,ii);
    % DONT ASSEMBLE... this will be faster
    %rhs1 = rhs - h*Md*v + h*Kd*v + h^2*K*v;
    % ... and while it might not be pretty, it is likewise faster to reuse
    % the common CBv expression...
    if ~integrator.useFullAinv 
        
        Afun = @(v) qspcgHelperFullK3D(h,M,D,K,v);
        fext = h*(D*v) + (h*h)*(K*v);
        rhs = rhs + fext;
    elseif integrator.useFullAinv
        % This is an option to compare/test the use of the full Ainv
        % vs the unassembled version. Both should give the same solution
        % cloth simulations need the full K to inject the bending energy
        % hessian
        h2 = h*h;
        h2K = h2 * K;
        hD = h* D;
        A = M - hD - h2K;     
        Afun = @(v) A*v;

        fext = hD*v + h2K*v;
        rhs = rhs + fext;
    end
    
    dofCount = 3 * sum([mesh3D.N]);
    
    cache.ApproximatedDeltaV = zeros(dofCount, 1);
    contactDv = zeros(size(rhs));

    % presolve on new contacts if we have any!
    if ~isempty( cache.cInfo ) 

        %force rigids in new contacts
%         rigidVerts = mesh3D.rigidIDbyVert(contactVertices) >= 1;
%         isNewContact(rigidVerts) = true;

        %force all contacts as new
%         isNewContact(:) = true;

        % Hmm... (sign may not make complete sense, but probably correct
        % based on what was done previously... note that we only use deltav
        % as the result of the PGS solve, so this is the main place where
        % lambdas are needed)
        tmpRhs = Jc' * warmStartLambdas;
        rhs = rhs + tmpRhs(ii); % existing forces on the RHS

        if any( isNewContact )
        
            % [  A  Jcn' ] [   dv   ]  = [    rhs   ]
            % [ Jcn  0   ] [ lambda ]    [  -Jcn v  ]

            % The key here is that we want Jcn v = 0, and solving for dv
            % this means Jcn dv = -Jcn v for the current velocity v.
            % So when we form the schur complement of this system we get
            % Jcn Ainv Jcn' lambda = Jcn Ainv rhs + Jcn v 
            % So that's what we need on the new contact rhs, ncRHS

            ncInds = find(isNewContact);
%             ncInds = 1:(size(Jc,1)/3);

            % Create a smaller Jc matrix for new contacts
            % Perhaps no need to solve for tangents?  only normals here!
           constraintsIDs = ncInds*3-2;
           Jcn = Jc( constraintsIDs, ii );
           %Todo: toggle untoggle tangents
%            JcnInds = reshape([ncInds*3-2,ncInds*3-1,ncInds*3]',1,[])';
%            Jcn = Jc( JcnInds, ii );


            if integrator.useFullContactAinv && integrator.useFullAinv %also needs full Ainv to be toggled on
                 % -----debugging version
                %uses the full A
                JcnAinv = Jcn/A;    
            else
                % We will use precomptued diagonal 3x3 blocks of A rather than 
                % the proper A inverse (or LDLT factorization)
                JcnAinv = Jcn * cache.AInvBlocks(ii,ii);
            end
            AtoSolve = JcnAinv * Jcn';
            contactVelocities = reshape([cInfo(ncInds).velocity],[],1);
            ncRHS = JcnAinv * rhs + Jcn * v + contactVelocities(1:3:end);%only normals in the contact frame
            ncRHS = ncRHS+integrator.Baumgarte * phi(ncInds);

            % small and dense (only new contacts) so just do a direct solve
            % or just let Matlab decide what's best with mldivide.
            if integrator.PGSquicksolve
                mu = [ cInfo(:).frictionCoefficient ]; %this is just to test, a real pgs implementation would skip the tangents parts.
                [ newLambdas ] = solveLDLTPGS3DwithJAinvJTHelper( settings.PGSiterations, zeros(numel(ncRHS),1), AtoSolve,ncRHS, mu, integrator.Compliance );
            else
                AtoSolve = AtoSolve - eye(size(AtoSolve))*integrator.Compliance;
                newLambdas = AtoSolve \ ncRHS;
%                 newLambdas(newLambdas<0) = 0;
            end
            % NOTE: h is baked into lambda!  See the equation above and 
            % note we are solving for dv, NOT acceleration, so lambda is
            % an impulse!
            rhs = rhs - Jcn' * newLambdas;

            if (integrator.useFullAinv || settings.addBendingEnergy || settings.useGrinspunPlanarEnergy) && (integrator.useQuicksolveContactFilter == 1||integrator.useQuicksolveContactFilter == 7)
                %fast neighbor filter contact dv version
                isContactless = true(mesh3D.N,1);
                isContactless(contactVertices) = false;
                newContactVertices = contactVertices(isNewContact);
                [newContactNeighbors,neighborLambdaID] = find(mesh3D.VertexAdjacencyMatrix(:,newContactVertices));
                addConstraint = isContactless(newContactNeighbors);
                normals = reshape([cInfo.normal],3,[])';
                newNormals = -normals(isNewContact,:);
                weightedLambdas = newLambdas(neighborLambdaID);

                %remove contactless constraints
                weightedLambdas = weightedLambdas(addConstraint);
                newContactNeighbors = newContactNeighbors(addConstraint);
                neighborLambdaID = neighborLambdaID(addConstraint);
                
                iiC = repmat(1:numel(weightedLambdas),3,1)';
                jjC = [3*newContactNeighbors-2,3*newContactNeighbors-1,newContactNeighbors*3];
                valsC = newNormals(neighborLambdaID,:);
                JcnNeighbor = sparse(iiC,jjC,valsC,numel(weightedLambdas),mesh3D.N*3);
                JcnNeighbor = JcnNeighbor(:,ii);
                JcnUp = [Jcn; JcnNeighbor];

                filterCenterWeight = 1;
                if integrator.useQuicksolveContactFilter == 7
                    filterCenterWeight = 4;
                end

                contactDv = JcnUp'*([filterCenterWeight*newLambdas;weightedLambdas]);
            elseif (integrator.useQuicksolveContactFilter == 3 || integrator.useQuicksolveContactFilter == 6)
                contactDv = A\(Jcn'*(newLambdas)); %actual contact dv(slow)
            elseif integrator.useQuicksolveContactFilter == 4
                contactDv = cache.AInvBlocks(ii,ii)*(Jcn'*(newLambdas)); %fast approx contact dv solve
            elseif integrator.useQuicksolveContactFilter == 5
                contactDv = cache.Lpre(ii,ii)*(Jcn'*(newLambdas));
            elseif integrator.useQuicksolveContactFilter == 8
                [contactDv,fl0,rr0,it0,rv0] = pcg(A,(Jcn'*(newLambdas)),1e-8,1,cache.Lpre',cache.Lpre);
            elseif integrator.useQuicksolveContactFilter == 10
                contactDv = cache.Lpre(ii,ii)\(Jcn'*(newLambdas));
            end
        end

        if settings.quicksolveSimulation %probably want to couple this with full new contacts when debugging
            cache.prevLambdas = warmStartLambdas;
        end
    end


    % final deltav is solved deltav from contacts + deltav from external forces.
    % This is always 1 iteration, the other version is also there to test
    % convergence
%     warmStartDv = zeros(size(rhs));
    warmStartDv = contactDv;
    
    if integrator.useQuicksolveContactFilter ==2  && any( isNewContact )
        warmStartDv = solveCGNoLoop(Afun,rhs,warmStartDv) ;
    end

    if settings.PCGiterations == 1 %This version should be much faster for 1 iteration by avoiding wasteful computations
        cache.ApproximatedDeltaV(ii) = PCGnoLoop(Afun, rhs, cache, mesh3D.objectDOFsPinned, warmStartDv);
%         [cache.ApproximatedDeltaV(ii),fl0,rr0,it0,rv0] = pcg(Afun, rhs,1e-10,1,cache.Lpre,cache.Lpre',warmStartDv);
    else
        [cache.ApproximatedDeltaV(ii),fl0,rr0,it0,rv0] = pcg(Afun, rhs,1e-10,settings.PCGiterations,cache.Lpre',cache.Lpre,warmStartDv);
    end
    cache.ApproximatedDeltaV(ii) = cache.ApproximatedDeltaV(ii) + gravityDv;

    if integrator.useQuicksolveContactFilter == 6
        cache.ApproximatedDeltaV(ii) = A\rhs + gravityDv;
    end

    %scripted animations positions
%     for i=1:numel(animationScripter)
%         [meshp, meshv] = animationScripter{i}.scriptPositions(mesh3D.p, mesh3D.v + cache.ApproximatedDeltaV, frame, h);
%         cache.ApproximatedDeltaV = meshv - mesh3D.v;
%     end

    %add animations
%     for i=1:numel(animationScripter)
%         cache.ApproximatedDeltaV = animationScripter{i}.scriptVelocityQuicksolve(mesh3D.p, mesh3D.v, cache.ApproximatedDeltaV, frame, h);
%     end

    if any(isnan(cache.ApproximatedDeltaV))
        errorQuicksolveNan = true
    end
    
    %debug option that allows to visualize the results of the quicksolve
    if settings.quicksolveSimulation %
        %this call is needed for the rigidificator to have an up to date
        %matrix when testing the quicksolve
        [~] = mesh3D.getB(cache);
        cache.oldp = mesh3D.p;
        cache.oldv = mesh3D.v;

        mesh3D.v = mesh3D.v + cache.ApproximatedDeltaV;
        mesh3D.v(mesh3D.pinnedDOFs) = mesh3D.v(mesh3D.pinnedDOFs) * 0;
        mesh3D.p = mesh3D.p + h * mesh3D.v;
        cache.clearWarmStartInfo();
    end

end