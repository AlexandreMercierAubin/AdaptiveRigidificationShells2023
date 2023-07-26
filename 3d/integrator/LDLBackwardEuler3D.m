classdef LDLBackwardEuler3D < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        % none needed
        projectToSPD = 0;
        prevAcontact
        prevA
        prevJc
        prevLambda
        recordConditionNumber = 0;
        conditionNumberList = [];
        conditionNumberContactList = [];
        recordPCGvsLDL = 0;
        PCGvsLDL = [];
    end
    
    methods
        function obj = LDLBackwardEuler3D()
            obj@Integrator();
            obj.Name = 'Backward Euler for 3D';
        end
        
        function integrate( obj, mesh3D, h, Jc, phi, cInfo, settings, cache, td, animationScripter, frame, energyModel,~)
            if nargin < 4
                Jc = zeros( 0, mesh3D.N*3 ); % no constraints
            end          
            if nargin < 5
                phi = [];
            end
            
            ticIntegrateForces = tic;
            
            [bigB, bigGamma] = mesh3D.getB(cache);
            Jc = Jc * bigGamma;
           
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
            
            %scripted animations force impulses
            for i=1:numel(animationScripter)
                mesh3D.f = animationScripter{i}.scriptForceAnimation(mesh3D.f, frame, h);
            end
            
            if isa( mesh3D, 'AdaptiveMesh3D' )
                mesh3D.computeRigidForce(h);
            end
            
            % compute right hand side
            rhs = h * mesh3D.getCurrentForce();
            if any(isnan(rhs))
                error =1
            end
            
            K = bigGamma'*cache.K*bigGamma;
            Kd = bigGamma'*cache.Kd*bigGamma;

            v = mesh3D.getCurrentVelocity();
            M = mesh3D.getM();
      
            D = bigGamma' * cache.D * bigGamma;
            hD = h*D;
            h2K=(h*h)*K;
            rhs = rhs + hD*v + h2K*v;
            A = M - hD - h2K;
            cache.A = A;
            
            if obj.projectToSPD
%                 A = sparse(nearestSPD(full(A)));
%                 A = nearestSPDS(A);
            end

            % make set of total active dofs
            ii = mesh3D.activeDOFs;
            tldl = tic;
            [L, D, P, S] = ldl(A(ii, ii));
                        
            % part of deltav that comes from outside forces
            % this is basically Ainv * rhs but without computing Ainv
            extDv = S * (P * (L' \ (D \ (L \ (P' * (S * rhs(ii)))))));
            timeldl = toc(tldl);
            if obj.recordPCGvsLDL
                G = bigGamma(mesh3D.unpinnedDOFs,ii);
                localLpre = tril(S*G'*cache.Lpre*G*S);
                pre = @(r) (localLpre'\(localLpre\(r)));
                Ascaled = S*A(ii,ii)*S;
                bscaled = S*rhs(ii);
%                 Afun = @(v) Ascaled*v;

                tpcg = tic;
                [extDv2,flag,res] = pcg( Ascaled, bscaled, 1e-4, 20, pre);
                timepcg = toc(tpcg);    
                
                obj.PCGvsLDL = [obj.PCGvsLDL;timeldl,timepcg];
            end

            if obj.recordConditionNumber == 1 %hummm use backwardEuler3D for this one without the scaling from LDL
                obj.conditionNumberList = [obj.conditionNumberList;cond(full(A))];
            elseif obj.recordConditionNumber == 3 %MC30 scaling
                obj.conditionNumberList = [obj.conditionNumberList;cond(full(L*D*L'))];
            end

            td.integrateForces = toc( ticIntegrateForces );
            ticIntegrateContact = tic;
            if isa(mesh3D, "AdaptiveMesh3D")
                numdofs = 3 * size(mesh3D.ElasticInds, 2) + 6*numel(mesh3D.RigidBodies);
                deltav = zeros( numdofs, 1);
            else
                deltav = zeros( mesh3D.N*3, 1);
            end
            
            if ( isempty(cInfo) || numel(ii)==0)
                deltav(ii) = extDv;
            else      
                Jca = Jc(:, ii);
                assert( issparse( Jca ) );
                Lrhs = Jca * (v(ii)+extDv); % lower right hand side           
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
                        warmStartLambdas = cache.findWarmStartLambdaArray();
                    end
                end
 
%                 Acontact = -Jc(:, ii)*(A(ii,ii)\Jc(:, ii)');
%                 if obj.projectToSPD
% %                 Acontact = sparse(nearestSPD(full(Acontact)));
% %                 Acontact = nearestSPDS(Acontact);
%                 end
%                 lambda = solvePGS3D(settings.PGSiterations, Acontact, Lrhs, warmStartLambdas, cache.cInfo);
%                 obj.prevAcontact = Acontact;
%                 obj.prevA = A;
%                 obj.prevJc = Jc;
%                 obj.prevLambda = lambda;
%                 dv = (A(ii,ii) \ (Jc(:, ii)' * lambda));
%                 [lambda, dv] = solveLDLTPGS3D(settings.PGSiterations, Jc(:, ii), L, D, P, S, Lrhs, warmStartLambdas, cInfo, obj.Compliance, td);
                % The with JAinvJT version is more than 2x faster... not
                % clear it ever makes sense to use the T version commented
                % above.
                
                JcaT = Jca';     % Sparse matrix transpose is not free!  probably cheap
                
                JcAinvJcT = Jca * (S * (P * (L' \ (D \ full(L \ (P' * (S * (JcaT))))))));
                cache.JcAinvJcT = JcAinvJcT;
%                 lambda = -JcAinvJcT\Lrhs;
                [lambda] = solveMexPGS3D(settings.PGSiterations, JcAinvJcT, Lrhs, warmStartLambdas, cInfo, obj.Compliance, td);
                dv =  S * (P * (L' \ (D \ full(L \ (P' * (S * (JcaT * lambda)))))));

                cache.prevLambdas = lambda;

                % final deltav is solved deltav from contacts + deltav from
                % external forces
                deltav(ii) = dv + extDv;

%                 if any(abs(lambda) > 1e6) %likely about to explode      
%                     figure
%                     imagesc(A-obj.prevA);
%                     colorbar;
%                     figure
%                     imagesc(Jc-obj.prevJc);
%                     colorbar;
%                     figure
%                     imagesc(Acontact-obj.prevAcontact);
%                     colorbar;
%                     printExplosion = 1
%                     
%                     indices = find(abs(lambda) > 1e4);
%                     contactIndices = floor((indices-1) ./ 3) + 1;
%                     pos = reshape([cInfo(contactIndices).point],3,[])';
%                     plotDebugSpheresAtPos(pos);
%                 end
            end

            if any(abs(deltav)>=1e5)
                warningHighVelocityChange = 1
%                 indices = find(abs(deltav)>1e5);
%                 acc = h*(mesh3D.M\mesh3D.f);
%                 deltav(indices) = acc(indices);
%                 settings.DrawVertexNumbers = unique(floor(indices./3))+1;
            end
            % store contact information for warm starts on next run
            cache.prevCInfo = cInfo;
            % the warm start cache has been used, can clear it now
            cache.clearWarmStartInfo;
            
            td.integrateContacts = toc( ticIntegrateContact );

            cache.oldp = mesh3D.p;
            cache.oldv = mesh3D.v;
            
            mesh3D.updateParticles(h, deltav(ii));

            %scripted animations positions
            for i=1:numel(animationScripter)
                [mesh3D.p, mesh3D.v] = animationScripter{i}.scriptPositions(mesh3D.p, mesh3D.v, frame, h);
            end
        end
    end
end