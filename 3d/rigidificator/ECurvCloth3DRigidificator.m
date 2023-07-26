classdef ECurvCloth3DRigidificator < Rigidificator
    %EDotMexRigidificator Uses EDot values to tell if a triangle should be
    %made rigid or not. This version uses mex to speedup the process
    % noting some interesting papers for this part
    % https://www.cs.columbia.edu/cg/pdfs/143-rods.pdf
    % https://www.cs.jhu.edu/~misha/Fall09/8-shells.pdf
    % https://studios.disneyresearch.com/wp-content/uploads/2019/03/Discrete-Bending-Forces-and-Their-Jacobians-Paper.pdf
    % https://gfx.cs.princeton.edu/pubs/Rusinkiewicz_2004_ECA/curvpaper.pdf
    % https://people.eecs.berkeley.edu/~jrs/meshpapers/MeyerDesbrunSchroderBarr.pdf
    properties
        ScaleByMaxEdgeLength
        RigidificationBendThreshold = 1e-3;
        ElastificationBendThreshold = 1e-2;
        bendType = 1; %1 for curvatures 2 for angles 3 for curvature from rest dual edge
        useAddNormalF = false;
        useCachedF = true;
        Fprev;
        recordCurvatures = false;
        curvatures = [];
    end
    
    methods
        function obj = ECurvCloth3DRigidificator()
            %EDOTVECTORRIGIDIFICATOR Same as EDotRigidficator but
            %vectorized and thus much faster
            obj@Rigidificator();
        end

        function setBendingThresholdsFromPlanar(obj,d)
            % when scale consistency is needed divide by the mesh diameter
            % d
            if nargin < 2
                d = 1;
            end
            obj.RigidificationBendThreshold = 4*sqrt(obj.RigidificationThreshold)/d;
            obj.ElastificationBendThreshold = 4*sqrt(obj.ElastificationThreshold)/d;
        end

        function [bendRate, currentBend, angles] = computeBendRates(obj, meshTri3D, p, prevCurvature, h)
            [currentBend, angles, ~] = computeCurvature(meshTri3D,p,1,obj.bendType);
            bendRate = abs(currentBend - prevCurvature)./h;
        end

        function [bendRate, currentBend, angles] = computeBendRatesFromAngles(obj, meshTri3D, p, prevCurvature, h, angles)
            [currentBend, ~] = computeCurvatureFromAngles(meshTri3D,p,1,obj.bendType,angles);
            bendRate = abs(currentBend - prevCurvature)./h;
        end
        
        function checkForElastification( obj, meshTri3D, cache, frame, h, settings )
            if ~isa(meshTri3D, 'AdaptiveMesh3D') 
                return;
            end
            
            frameIndex = mod(frame, obj.FrameCount) + 1;
            %rigidification Edot============================
            % init the Edots
            if ~numel(meshTri3D.RigidificationValues)
                  meshTri3D.RigidificationValues = zeros(size(meshTri3D.t, 1), obj.FrameCount);
%                   obj.Fprev = cache.F;
                  obj.Fprev = mexComputeFtri3D(meshTri3D.getPositionFormatted,meshTri3D.t,meshTri3D.triNormals,meshTri3D.dphidx(:));
            end
            if ~obj.useAddNormalF && ~obj.useCachedF
                F = mexComputeFtri3D(meshTri3D.getPositionFormatted,meshTri3D.t,meshTri3D.triNormals,meshTri3D.dphidx(:));
    %             F = computeFtri3D(meshTri3D.getPositionFormatted,meshTri3D.t,meshTri3D.triNormals,meshTri3D.dphidx);
            else
                F = cache.F;
            end
            
            %computing this here because we need to update the prev angles anyway
            if settings.addBendingEnergy % use the angles computed from the bending energy to save on some computations
                [curvatureRate, currentCurvature, angles] = obj.computeBendRatesFromAngles(meshTri3D, meshTri3D.p, meshTri3D.prevCurvature,h, cache.dihedralAngles);
            else
                [curvatureRate, currentCurvature, angles] = obj.computeBendRates(meshTri3D, meshTri3D.p, meshTri3D.prevCurvature,h);
            end
            
            if obj.recordCurvatures
                obj.curvatures = [obj.curvatures,currentCurvature];
            end

            cache.bendingRate = curvatureRate;
            if isempty(meshTri3D.ElasticTetInds) 
                meshTri3D.RigidificationValues(:, frameIndex) = zeros(size(meshTri3D.t, 1),1);
                rigidTrisHist = meshTri3D.RigidificationValues < obj.RigidificationThreshold;
                cache.edotnorms = zeros(size(meshTri3D.t,1),1);
            elseif ~settings.RigidificationEnabled || ~meshTri3D.EnableAutomaticRigidification
                mesh2d.RigidificationValues(:, frameIndex) = ones(size(meshTri3D.t, 1),1)*obj.RigidificationThreshold;
                rigidTrisHist = mesh2d.RigidificationValues < obj.RigidificationThreshold;
                cache.edotnorms = [];
            else
                EDotNorms = mexEdiffNorm3D( obj.Fprev, F, h);
%                 EDotNorms = EdiffNorm( F, obj.Fprev, h,3);
                cache.edotnorms = EDotNorms; % for plotting
                
                % set EDots
                logicalIndices = false(size(meshTri3D.t,1),1);
                logicalIndices(meshTri3D.ElasticTetInds) = true;
                meshTri3D.RigidificationValues(logicalIndices, frameIndex) = EDotNorms(logicalIndices);%TODO: make sure this doesn't cause an indexing issue
                meshTri3D.RigidificationValues(~logicalIndices, :) = 0;
                %for tri3D meshes, look at the bending
                cache.bendingRate = curvatureRate;
                %TODO: the curvature rates are big for the disk
                triBend = accumarray([meshTri3D.bendingEdgesElements(:,1);meshTri3D.bendingEdgesElements(:,2)],[curvatureRate;curvatureRate],[size(meshTri3D.t,1),1], @max);

                edgeBendTest = all(triBend < obj.RigidificationBendThreshold,2);
                meshTri3D.RigidificationBendValues(:,frameIndex) = edgeBendTest;
                rigidBend = all(meshTri3D.RigidificationBendValues,2);

                rigidStretch = all(meshTri3D.RigidificationValues < obj.RigidificationThreshold,2);
                rigidTrisHist = rigidBend&rigidStretch;

                if ~settings.RigidificationEnabled
                    rigidTrisHist(:,frameIndex) = false(size(meshTri3D.RigidificationValues));
                end
            end
            rigidTrisHist(meshTri3D.AlwaysRigidElements,:) = true;
            %End rigidification
            %============================================
            %Elastification Edot ====================
            if isempty(meshTri3D.RigidBodies) || ~settings.ElastificationEnabled || ~meshTri3D.EnableAutomaticRigidification
                trisToElastify = false(size(meshTri3D.t,1),1);
            else
                rowsIdRange = 1:size(meshTri3D.t, 1) * 9;
                setDiffIds = true(numel(rowsIdRange),1);
                setDiffIds(meshTri3D.ActiveBRows) = false;
                rigidRows = rowsIdRange(setDiffIds);

                trisIdRange = 1:size(meshTri3D.t, 1);
                setDiffIds = setDiffIds(1:9:end);
                rigidTris = trisIdRange(setDiffIds);

                vel = meshTri3D.v + cache.ApproximatedDeltaV; 

                approxPositions = meshTri3D.p + h*vel;
                approxNormals = mexNormals(meshTri3D.formatPositions(approxPositions),meshTri3D.t);
                if obj.useAddNormalF
                    [steppedB, steppedF]=addShellNormalDeformation(meshTri3D, approxPositions, approxNormals, settings); % not sure if this is required
                else 
                    steppedF = mexComputeFtri3D(meshTri3D.formatPositions(approxPositions),meshTri3D.t,approxNormals,meshTri3D.dphidx(:));
%                 steppedF = computeFtri3D(meshTri3D.formatPositions(approxPositions),meshTri3D.t,approxNormals,meshTri3D.dphidx);
                end
                EDotApproxNorms = mexEdiffNorm3D( F(rigidRows), steppedF(rigidRows),h );
%                 EDotApproxNorms = EdiffNorm( F(rigidRows), steppedF(rigidRows),h, 3);
               
                cache.edotApproxNorms = EDotApproxNorms; % for plotting

                % rigid tris with value > threshold in Rigidification causes
                % issues, so this prevents that. This might not be necessary
                % though.
                meshTri3D.RigidificationValues(rigidTris, :) = 0;
                
                % edot norms again but back into a vector of the entire mesh
                % triangles rather than just the rigids
                globalEDots = zeros(size(meshTri3D.t, 1),1);
                globalEDots( rigidTris ) = EDotApproxNorms;

                % new elastic triangles
                elasticStretch = globalEDots > obj.ElastificationThreshold;

                %bending angle test
                [curvatureRateApprox,~] = obj.computeBendRates(meshTri3D, approxPositions, currentCurvature,h);
                cache.dihedralApproxRate = curvatureRateApprox;
                triBend = accumarray([meshTri3D.bendingEdgesElements(:,1);meshTri3D.bendingEdgesElements(:,2)],[curvatureRateApprox;curvatureRateApprox],[size(meshTri3D.t,1),1], @max);
                triBend(meshTri3D.isTetElastic) = 0;
                edgeBendTest = any(triBend > obj.ElastificationBendThreshold,2);

                trisToElastify = elasticStretch|edgeBendTest;
                
                %This is NOT needed as the next step rigidification check
                %will have high strain rates if needed because it will be
                %simulated as elastic due to trisToElastify being toggled
%                 meshTri3D.RigidificationValues(trisToElastify, frameIndex) = globalEDots(trisToElastify);
                
            end
            trisToElastify(meshTri3D.AlwaysRigidElements) = false;
            %End elastification
            %============================================
            %update prevdihedralangle
            cache.prevPrevDihedralAngle = angles;
            meshTri3D.prevDihedralAngle = angles;
            meshTri3D.prevCurvature = currentCurvature;
            obj.Fprev = F;

            numRigid = 0;
            if frame >= obj.FrameCount || settings.FirstFrameRigidification
                if obj.PreventPinnedRigidification
                    trisToElastify(meshTri3D.pinnedTets) = true;
                end
                if ~isempty(meshTri3D.animationInds)   
                    animated = false(numel(meshTri3D.pinned),1);
                    animated(meshTri3D.animationInds) = true;
                    pinned = meshTri3D.pinned;
                    pinned(meshTri3D.animationInds) = true;
                    pinned = logical(pinned);
                    pinnedVertPerTri = pinned(meshTri3D.t(:,1:3));
                    stablePinnedElements = find(sum(pinnedVertPerTri,2) >= 2);
                    isStableTri = false(size(meshTri3D.t,1),1);
                    isStableTri(stablePinnedElements) = true;
                else
                    pinned = meshTri3D.pinned;
                    pinned = logical(pinned);
                    stablePinnedElements = meshTri3D.stablePinnedElements;
                    isStableTri = false(size(meshTri3D.t,1),1);
                    isStableTri(stablePinnedElements) = true;
                    animated = false(numel(meshTri3D.pinned),1);
                end
                
                [numRigid, rigidIDbyElement, rigidIDbyVert, isVertElastic, isVertBoundary, isElementElastic, isComponentPinned, isAnimatedComponent] = mexRigidBodyConnectedMixed(rigidTrisHist, trisToElastify, meshTri3D.AdjagencyMatrix, meshTri3D.t, meshTri3D.N , pinned, meshTri3D.valence, stablePinnedElements, isStableTri, animated);
                meshTri3D.rigidIDbyVert = rigidIDbyVert;

%                 stateChanged = isElementElastic~=meshTri3D.isTetElastic;
%                 meshTri3D.RigidificationValues(stateChanged & ~meshTri3D.isTetElastic, :) = inf;
                
                hasNotChanged = all(isElementElastic==meshTri3D.isTetElastic);
                if hasNotChanged
                    return;
                end
                meshTri3D.isTetElastic = isElementElastic;
                
                %disgusting
                meshTri3D.ElasticTetInds = find(isElementElastic)';
                meshTri3D.ElasticInds = find(isVertElastic)';
                
                rigidVert = find(~isVertElastic);
                inds2 = rigidVert * 2;
                RigidDOFs = reshape([inds2 - 1; inds2], 1, []);
                
                rigidVert = find(~isVertElastic);
                
                inds3 = meshTri3D.ElasticInds * 3;
                inds9 = meshTri3D.ElasticTetInds * 9;
                meshTri3D.ElasticDOFs = reshape([inds3 - 2; inds3 - 1; inds3], 1, []);
                meshTri3D.ActiveBRows = reshape([inds9 - 8; inds9 - 7; inds9 - 6; inds9 - 5; inds9 - 4; inds9 - 3; inds9 - 2; inds9 - 1; inds9], 1, []);
    
                [ com, comdot, mass, rotMass, angularMomentum, vertexDisp] = mexRigidBodyProperties3D( numRigid, rigidIDbyVert, meshTri3D.p, meshTri3D.v, meshTri3D.mass );

                inds3 = rigidVert * 3;
                meshTri3D.VertexDisp(rigidVert,:) = [vertexDisp(inds3-2), vertexDisp(inds3-1), vertexDisp(inds3)];
            end
            
            meshTri3D.RigidBodies(numRigid + 1:end) = [];
            lenRigid = numel(meshTri3D.RigidBodies);
            for i = 1:numRigid
                if lenRigid < i
                    meshTri3D.RigidBodies(i) = RigidBody3D(meshTri3D);
                end
                meshTri3D.RigidBodies(i).TetInds = find(rigidIDbyElement == i)';
                meshTri3D.RigidBodies(i).VertexIDs = find(rigidIDbyVert == i)';
                meshTri3D.RigidBodies(i).Position = com(:,i);
                meshTri3D.RigidBodies(i).Velocity = comdot(:,i);
                meshTri3D.RigidBodies(i).Mass = mass(i);
                inertia = rotMass(:,i);
                meshTri3D.RigidBodies(i).Inertia = reshape(inertia,3,3);
                meshTri3D.RigidBodies(i).Inertia0 = meshTri3D.RigidBodies(i).Inertia;
                meshTri3D.RigidBodies(i).AngularVelocity = meshTri3D.RigidBodies(i).Inertia \ angularMomentum(:,i);
                meshTri3D.RigidBodies(i).Rotation = eye(3);
                
                inds3= meshTri3D.RigidBodies(i).VertexIDs * 3;
                meshTri3D.RigidBodies(i).DOFs = reshape([inds3 - 2; inds3 - 1; inds3], 1, []);
                meshTri3D.RigidBodies(i).Force = [0;0;0];
                meshTri3D.RigidBodies(i).Torque = 0;
                meshTri3D.RigidBodies(i).isPinned = isComponentPinned(i);
                meshTri3D.RigidBodies(i).isAnimated = false; %Todo: implement scripted animations for cloth
            end

            if frame >= obj.FrameCount || settings.FirstFrameRigidification
                elasticDOFsCount = numel(meshTri3D.ElasticDOFs);
                n = elasticDOFsCount + numel(meshTri3D.RigidBodies) * 6;
                extendedMdiag = zeros(n, 1);
                extendedMdiag(1:elasticDOFsCount) = meshTri3D.mass(meshTri3D.ElasticDOFs);

                if ~isempty(meshTri3D.RigidBodies)
                    extendedMdiag(elasticDOFsCount + 1:6:end) = [meshTri3D.RigidBodies.Mass];
                    extendedMdiag(elasticDOFsCount + 2:6:end) = [meshTri3D.RigidBodies.Mass];
                    extendedMdiag(elasticDOFsCount + 3:6:end) = [meshTri3D.RigidBodies.Mass];
                end

                meshTri3D.AdaptiveM = spdiags(extendedMdiag, 0, n, n);

                for i = 1:numel(meshTri3D.RigidBodies)
                    inds = elasticDOFsCount + (i-1)*6 + (4:6);
                    meshTri3D.AdaptiveM(inds,inds) =  meshTri3D.RigidBodies(i).Inertia;
                end
                
                meshTri3D.computeActiveDOFs();
            end
        end
    end
end

