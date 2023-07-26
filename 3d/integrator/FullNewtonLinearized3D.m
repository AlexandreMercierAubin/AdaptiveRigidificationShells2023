classdef FullNewtonLinearized3D < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        % none needed
        prevWarmstarts
        regularizator = 1;
        useGamma = 0;
        dphidX
        solverType
        solverTypeEnum
        projectToSPD = 0;
        useLineSearch = true;
        maxIterations = 5;
        eps = 1e-9;
    end
    
    methods
        function obj = FullNewtonLinearized3D()
            obj@Integrator();
            obj.Name = 'Newton solve Backward Euler for 3D';
            obj.maxIterations = 5;
        end
        
        function integrate( obj, mesh3D, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame, energyModel, collisionDectectors )
            if nargin < 4
                Jc = sparse( 0, mesh3D.N*3 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            
            if size(mesh3D.t,2) == 3 && settings.addShellNormalDeformation
                dphidX = reshape(mesh3D.dphidx(:)',9,size(mesh3D.t,1))';
            else
%                 dphidX = reshape(mesh3D.dphidx(:)',12,size(mesh3D.t,1))';
%                 dphidX =  linear_tetmesh_dphi_dX(mesh3D.getPosition0Formatted,double(mesh3D.t));
            end
            ticIntegrateForces = tic;
            ii = mesh3D.activeDOFs;
            
            [bigB, bigGamma] = mesh3D.getB(cache);
            JcAdaptive = sparse(Jc * bigGamma);
            
            mesh3D.resetForce;
            mesh3D.f = mesh3D.f + obj.CustomForce;
            for i=1:numel(animationScripts)
                mesh3D.f = animationScripts{i}.scriptForceAnimation(mesh3D.f, frame, h);
            end
            mesh3D.applyAcceleration([0; 0; obj.Gravity]);
            fext = mesh3D.f;
            
            if isa( mesh3D, 'AdaptiveMesh3D' )
                mesh3D.computeRigidForce(h);
                tetsInds = 1:size(mesh3D.t,1);
                rigidTets = tetsInds(mesh3D.isTetElastic)';
                rigidFentries = [rigidTets*9-8,rigidTets*9-7,rigidTets*9-6,rigidTets*9-5,rigidTets*9-4,rigidTets*9-3,rigidTets*9-2,rigidTets*9-1,rigidTets*9];
                rigidFinds = reshape(rigidFentries',1,[]);
            else
                rigidDofs = [];
            end
            
            v0 = mesh3D.getCurrentVelocity();%initial iterate
            v0ii = v0(ii);

            % compute right hand side
            M = mesh3D.getM();
            Mii = M(ii,ii);
            Minv = mesh3D.Minv;
            Minvii = Minv(ii,ii);
            Md = bigGamma' * mesh3D.Md * bigGamma;
            Mdii = Md(ii,ii);
            deltavfull = v0*0;

            if numel(deltavfull(ii)) > 0
                [deltavii, objective] = handmadeNewton(@energyMomentumLinearized,deltavfull(ii));
            end

            cache.oldp = mesh3D.p;
            cache.oldv = mesh3D.v;
            
%             deltavii = (vii-v0ii);
            mesh3D.updateParticles(h, deltavii);
            cache.oldDv = mesh3D.v - cache.oldv;

            %scripted animations positions
            for i=1:numel(animationScripts)
                [mesh3D.p, mesh3D.v] = animationScripts{i}.scriptPositions(mesh3D.p, mesh3D.v, frame, h);
            end

            update = toc(ticIntegrateForces);
            td.integrateForces = td.integrateForces + update;

            function [deltav, objective] = handmadeNewton(energyMomentum, deltavii);
                deltav = deltavii;
                elasticDv = deltav;
                sigma = 1e-4;
                giveUpLS = 1e-15;
                
                stop = false;
                iterations=0;
                while iterations < obj.maxIterations && ~stop
                    [objective,g,H] = energyMomentum(elasticDv);
                    if obj.projectToSPD
                        H = nearestSPDS(H);
                    end
                    cache.A = H;
                    newDeltav = -(H\g);
                    deltavUpdate = newDeltav-elasticDv;
                    alpha = 1;
                    if obj.useLineSearch
                        [steppedObjective] = energyMomentum(elasticDv + alpha*deltavUpdate);
                        while objective - steppedObjective < sigma * alpha * objective && steppedObjective >  obj.eps && alpha > giveUpLS
                            alpha = alpha/2;
                            [steppedObjective] = energyMomentum(elasticDv + alpha*deltavUpdate);
                        end
                    end
                    
                    elasticDv = elasticDv + alpha*deltavUpdate;

                    iterations = iterations +1;

                    if steppedObjective < obj.eps
                        break;
                    end
                end

                [Jc,phi, cInfos, bigGamma] = detectContacts(collisionDectectors, cache, elasticDv, frame + 1);
                contactDv = solveContacts(H,elasticDv, Jc, phi, bigGamma);
                deltav = contactDv + elasticDv;
            end

            function [psi] = tmpGrinspun(obj,mesh3D,cache, p)
                isShell = mesh3D.elementType == mesh3D.elementTypeEnum.Shell;
                [ iiG, jjG, vals, dpsidX , valsD, Wa, Wl] = mexGrinspunPlanar(p,mesh3D.t,mesh3D.elkl,mesh3D.elka,mesh3D.area,isShell, mesh3D.elAlpha1, mesh3D.grinspunEnergyEdgeRestLength, sum(isShell));
                cache.grinspunPlanarH = -sparse(iiG,jjG,vals,mesh3D.N*3,mesh3D.N*3);
                cache.grinspunPlanarHd = -sparse(iiG,jjG,valsD,mesh3D.N*3,mesh3D.N*3);
                cache.grinspunPlanarForces = -dpsidX;
                allpsi = Wa + Wl;
                sumWa = sum(Wa);
                sumWl = sum(Wl);
                psi = sumWa + sumWl;
            end

            function [psi] = computeElasticEnergy(obj,mesh3D,cache, p)
                if settings.useGrinspunPlanarEnergy
                    isShell = mesh3D.elementType == mesh3D.elementTypeEnum.Shell;
                    ka = repmat(settings.ka,size(mesh3D.t,1),1);
                    kl = repmat(settings.kl,size(mesh3D.t,1),1);
                    [ iiG, jjG, vals, dpsidX , valsD, Wa, Wl] = mexGrinspunPlanar(p,mesh3D.t,kl,ka,mesh3D.area,isShell, mesh3D.elAlpha1, mesh3D.grinspunEnergyEdgeRestLength, sum(isShell));
                    cache.grinspunPlanarH = -sparse(iiG,jjG,vals,mesh3D.N*3,mesh3D.N*3);
                    cache.grinspunPlanarHd = -sparse(iiG,jjG,valsD,mesh3D.N*3,mesh3D.N*3);
                    cache.grinspunPlanarForces = -dpsidX;
                    allpsi = Wa + Wl;
                    sumWa = sum(Wa);
                    sumWl = sum(Wl);
                    psi = sumWa + sumWl;
                else
                    energyModel.computeEnergy(mesh3D,F);
                    psi = sum(energyModel.energy);
                end
            end

            function [objective,g,H]=energyMomentumLinearized(deltav) 
                
                v = v0;
                v(ii) = v0ii+deltav;
                vii = v(ii);
                [steppedPositions, angularVelocity, rigidBodyPosition, inertia] = mesh3D.peekNextParticles(h, deltav);
%                 steppedPositions= mesh3D.p;
                
                steppedPositions = SVDStrainLimiting(mesh3D, steppedPositions, cache, settings);
                [mesh3D.B, cache.F]=addShellNormalDeformation(mesh3D, steppedPositions, settings);
                
                F = cache.F;

                if nargout > 1
                    elasticEnergy = computeElasticEnergy(obj,mesh3D,cache, steppedPositions);
                    [psiBending, cache.bendingForces,cache.shellBendingH] = computeBendingGradHess(mesh3D, steppedPositions, cache, settings);
                else
                    elasticEnergy = computeElasticEnergy(obj,mesh3D,cache, steppedPositions);
                    psiBending = computeBendingGradHess(mesh3D, steppedPositions, cache, settings);
                end

                momentumPreservation = 0.5*vii'*Mii*vii - vii'*Mii*v0ii;
                objective = momentumPreservation + elasticEnergy + psiBending - h*vii'*fext(ii);
                
                if nargout > 1
                    %need to update bigB as we changed B
                    [bigB, bigGamma] = mesh3D.peekB(cache, rigidBodyPosition', steppedPositions);
                    elasticForces = cache.elasticForces; %TODO figure out why cahce.elasticForces work but not the energyModel.elasticForces
                    bendingForces = cache.bendingForces;
                    
                    sumForces = elasticForces + bendingForces + fext;
                    forces = mesh3D.peekCurrentForce(h, angularVelocity, rigidBodyPosition, inertia, sumForces, steppedPositions);

                    C = energyModel.derivative2HessianC(mesh3D.ActiveBRows, mesh3D.ActiveBRows);
                        
                    K = bigB'*C*bigB;
                    bigAlpha1 = mesh3D.bigAlpha1( mesh3D.ActiveBRows, mesh3D.ActiveBRows );
                    Kd = bigB'*(bigAlpha1*C)*bigB;

                    bigH = bigGamma' * cache.shellBendingH * bigGamma;
                    K = K + bigH;

                    g =  -(h*forces(ii) - h*(Mdii*vii - Kd(ii,ii)*vii) + h^2*K(ii,ii)*vii); 

                    if nargout > 2
                        H = Mii - h * (-Mdii + Kd(ii,ii)) - h^2 * K(ii,ii);
                    end
                end
            end

            function [Jc,phi, cInfos, bigGamma] = detectContacts(collisionDetector, cache, deltav, frame)
%                 [p, angularVelocity, rigidBodyPosition, inertia] = mesh3D.peekNextParticles(h, deltav);
                %collision detection synergizes better with rigidification
                %if we solve it for the current positions and latest
                %velocity update
                [steppedPositions, angularVelocity, rigidBodyPosition, inertia] = mesh3D.peekNextParticles(h, deltav(ii));
                [bigB, bigGamma] = mesh3D.peekB(cache, rigidBodyPosition', steppedPositions);
                Jc = sparse( 0, mesh3D.N*3 );
                phi = [];
                cInfos = contactInfo3D.empty;
                for j = 1:numel(collisionDetector)
                    [ Jci, phii, cInfosi ] = collisionDetector{j}.findContacts(mesh3D, frame, steppedPositions);
                    Jc = [ Jc; Jci ];
                    phi = [ phi; phii];
                    cInfos = [ cInfos, cInfosi ];
                end
                cache.contactIDs = [ cInfos(:).contactID ];
                cache.cInfo = cInfos;
            end

            function contactDv = solveContacts(Aii, extDeltav, Jc, phi, bigGamma)
                % store contact information for warm starts on next run
                cache.prevCInfo = cInfo;
                % the warm start cache has been used, can clear it now
                cache.clearWarmStartInfo;

                if size(Jc,1) == 0
                    contactDv=extDeltav*0;
                    return
                end
                JcAdaptive = sparse(Jc * bigGamma);
                Jca = JcAdaptive(:, ii);
                Lrhs = Jca * (v0ii + extDeltav); % lower right hand side         
                %adding baumgarte to normal velocities
                Lrhs(1:3:end) = Lrhs(1:3:end) + obj.Baumgarte * phi;
                %but not the the friction
                % LRHS gets velocity term for contacts with moving obstacles

                n = size(Lrhs, 1);
                warmStartLambdas = zeros(n,1);
                if ( settings.WarmStartEnabled )
                    if ~isempty( cache.prevCInfo )
                        [warmStartLambdas,newContacts] = cache.findWarmStartLambdaArray();
                    end
                end
 
                Acontact = JcAdaptive(:, ii)*(Aii\JcAdaptive(:, ii)');
                
                obj.prevWarmstarts = warmStartLambdas;
                
                assert( issparse( JcAdaptive ) );
                cache.JcAinvJcT = Acontact;
                lambda = solveMexPGS3D(settings.PGSiterations, Acontact, Lrhs, warmStartLambdas, cache.cInfo, obj.Compliance, td);

                % final deltav is solved deltav from contacts + deltav from
                % external forces
                cache.prevLambdas = lambda;
                contactDv = v0*0;
                contactDv(ii) = (Aii \ (JcAdaptive(:, ii)' * lambda));

                % store contact information for warm starts on next run
                cache.prevCInfo = cInfo;
                % the warm start cache has been used, can clear it now
                cache.clearWarmStartInfo;
            end




        end
    end
end