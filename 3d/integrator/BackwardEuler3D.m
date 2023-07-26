classdef BackwardEuler3D < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        % none needed
        prevWarmstarts
        regularizator = 1;
        useGamma = 0;
        projectToSPD = 0;
        recordConditionNumber = 0;
        conditionNumberList = [];
        conditionNumberContactList = [];
    end
    
    methods
        function obj = BackwardEuler3D()
            obj@Integrator();
            obj.Name = 'Backward Euler for 3D';
        end
        
        function integrate( obj, mesh3D, h, Jc, phi, cInfo, settings, cache, td, animationScripts, frame, energyModel, ~)
            if nargin < 4
                Jc = zeros( 0, mesh3D.N*3 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            
            ticIntegrateForces = tic;
            ii = mesh3D.activeDOFs;
            
            [bigB, bigGamma] = mesh3D.getB(cache);
            JcAdaptive = sparse(Jc * bigGamma);
            
            mesh3D.resetForce;
            
            if settings.useGrinspunPlanarEnergy %assumes the elastic forces is not containing planar forces, only tets
                internalf = cache.grinspunPlanarForces;
            else
                internalf = cache.elasticForces;
            end

             % Gravity (z is vertical in matlab)
            mesh3D.f = zeros( mesh3D.N*3, 1 );
            mesh3D.applyAcceleration([0; 0; obj.Gravity]);
            mesh3D.f = mesh3D.f + internalf + obj.CustomForce;

            if (settings.addBendingEnergy)
                mesh3D.f = mesh3D.f + cache.bendingForces;
            end
            
            for i=1:numel(animationScripts)
                mesh3D.f = animationScripts{i}.scriptForceAnimation(mesh3D.f, frame, h);
            end
            
            if isa( mesh3D, 'AdaptiveMesh3D' )
                mesh3D.computeRigidForce(h);
            end
            
            % compute right hand side
            rhs = h * mesh3D.getCurrentForce;

            bigAlpha1 = mesh3D.bigAlpha1( mesh3D.ActiveBRows, mesh3D.ActiveBRows );
            
            % assembly is somewhat costly... perhaps can be a bit more
            % efficient about it...
            % 1) sparsity structure is fixed... this can be exploited
            % 2) K and Kd do not need to be built separately!
            
            K = bigGamma'*cache.K*bigGamma;
            Kd = bigGamma'*cache.Kd*bigGamma;
            
            v = mesh3D.getCurrentVelocity();
            
            Md = bigGamma' * mesh3D.Md * bigGamma;
            M = mesh3D.getM();
      
            rhs = rhs - h*Md*v + h*Kd*v + h^2*K*v;
            A = M - h * (-Md + Kd) - h^2 * K;

            % make set of total active dofs
            
            Aii = A(ii,ii);
            extDeltav = Aii \ rhs(ii);

            if obj.recordConditionNumber == 1
                obj.conditionNumberList = [obj.conditionNumberList;cond(full(A))];
            elseif obj.recordConditionNumber == 2
                diagA = diag(Aii);
                S = diag(1./diagA);
                obj.conditionNumberList = [obj.conditionNumberList;cond(full(S*Aii))];
            end
                        
            % part of deltav that comes from outside forces
            % this is basically Ainv * rhs but without computing Ainv
            
            td.integrateForces = toc( ticIntegrateForces );
            ticIntegrateContact = tic;
            if isa(mesh3D, "AdaptiveMesh3D")
                deltav = zeros( 3 * size(mesh3D.ElasticInds, 2) + 6*numel(mesh3D.RigidBodies), 1);
            else
                deltav = zeros( mesh3D.N*3, 1);
            end
            
            if ( isempty(cInfo) )
                deltav(ii) = extDeltav;
            else            
                Lrhs = JcAdaptive(:, ii) * (v(ii) + extDeltav); % lower right hand side         
                %adding baumgarte to normal velocities
                Lrhs(1:3:end) = Lrhs(1:3:end) + obj.Baumgarte * phi;
                %but not the the friction
                % LRHS gets velocity term for contacts with moving obstacles
                contactVelocities = reshape([cInfo.velocity],[],1);
                Lrhs = Lrhs + contactVelocities;

                n = size(Lrhs, 1);
                warmStartLambdas = zeros(n,1);
                if ( settings.WarmStartEnabled )
                    if ~isempty( cache.prevCInfo )
                        [warmStartLambdas,newContacts] = cache.findWarmStartLambdaArray();
                    end
                end
 
%                 Acontact = -JcAdaptive(:, ii)*Ainv*JcAdaptive(:, ii)';
                Acontact = -JcAdaptive(:, ii)*(Aii\JcAdaptive(:, ii)');
                if obj.regularizator == 1
                    compliance = eye(size(Jc,2));
                    Winv = Jc*compliance*Jc';
                elseif obj.regularizator == 2
                    W = mesh3D.Mii;
                    Winv = Jc*(W\Jc');
                elseif obj.regularizator == 3
                    W = cache.Lpre*cache.Lpre';
                    Winv = Jc*(W\Jc');
                elseif obj.regularizator == 4
                    W = cache.AInvBlocks;
                    Winv = Jc*(W*Jc');
                elseif obj.regularizator == 5
                    W = cache.Apre;
                    Winv = Jc*(W\Jc');
                elseif obj.regularizator == 6
                    elC = cache.C;
                    elAlpha1 = mesh3D.bigAlpha1;
                    elK = sparse(mesh3D.B' * elC * mesh3D.B);
                    elKd = sparse(mesh3D.B' * (elAlpha1 * elC) * mesh3D.B);
                    elM = mesh3D.M;

                    dampingTerm = h*(-mesh3D.Md + elKd);
                    W = elM  - h^2 * elK - dampingTerm;
                    Winv = Jc*(W\Jc');
                end
                
                if obj.useGamma == 1
                    gamma = obj.Compliance/eigs(Winv,1);
                elseif obj.useGamma == 0
                    gamma = 1;
                else
                    gamma = obj.useGamma;
                end

                Acontact = Acontact - gamma*Winv;
                if obj.recordConditionNumber
                    obj.conditionNumberContactList = [obj.conditionNumberContactList;cond(Acontact)];
                end
                
                obj.prevWarmstarts = warmStartLambdas;
                
                assert( issparse( Jc ) );
                lambda = solvePGS3D(settings.PGSiterations, Acontact, Lrhs, warmStartLambdas, cache.cInfo);

                cache.prevLambdas = lambda;

                % final deltav is solved deltav from contacts + deltav from
                % external forces
                cache.prevLambdas = lambda;
%                 deltav(ii) = Ainv *(JcAdaptive(:, ii)' * lambda) + extDeltav;
                deltav(ii) = (Aii \ (JcAdaptive(:, ii)' * lambda)) + extDeltav;
            end

            % store contact information for warm starts on next run
            cache.prevCInfo = cInfo;
            % the warm start cache has been used, can clear it now
            cache.clearWarmStartInfo;
            
            td.integrateContacts = toc( ticIntegrateContact );
            tic
            
            cache.oldp = mesh3D.p;
            cache.oldv = mesh3D.v;
            
%             for i = 1: numel(animationScripts)
%                 deltav = animationScripts{i}.scriptVelocity(mesh.p, mesh.v, deltav, mesh.ElasticDOFs, mesh.ActiveDofsCorrespondingID, frame, mesh.rigidIDbyVert, h);
%             end
            deltav = deltav(ii);
            
            mesh3D.updateParticles(h, deltav);
            cache.oldDv = mesh3D.v - cache.oldv;
            update = toc;
            td.integrateForces = td.integrateForces + update;
            
            %scripted animations positions
            for i=1:numel(animationScripts)
                [mesh3D.p, mesh3D.v] = animationScripts{i}.scriptPositions(mesh3D.p, mesh3D.v, frame, h);
            end

        end
    end
end