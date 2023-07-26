classdef symplectic3D < Integrator
    % CHOLESKIBACKWARDEULER Backward euler integrator that uses choleski to
    % solve for M - h^2K

    properties
        % none needed
        projectToSPD = 0;
        prevAcontact
        prevA
        prevJc
        prevLambda
    end
    
    methods
        function obj = symplectic3D()
            obj@Integrator();
            obj.Name = 'simplectic for 3D';
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
            
            acceleration = mesh3D.M\mesh3D.f;
            deltav = acceleration * h;

            if any(abs(deltav)>=10)
                errExplosion = 1
                indices = find(abs(deltav)>1e5);
%                 settings.DrawVertexNumbers = unique(floor(indices./3))+1;
            end
            % store contact information for warm starts on next run
            cache.prevCInfo = cInfo;
            % the warm start cache has been used, can clear it now
            cache.clearWarmStartInfo;
            
%             td.integrateContacts = toc( ticIntegrateContact );

            cache.oldp = mesh3D.p;
            cache.oldv = mesh3D.v;
            ii = mesh3D.activeDOFs;
            mesh3D.updateParticles(h, deltav(ii));

            %scripted animations positions
            for i=1:numel(animationScripter)
                [mesh3D.p, mesh3D.v] = animationScripter{i}.scriptPositions(mesh3D.p, mesh3D.v, frame, h);
            end
        end
    end
end