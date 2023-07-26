classdef Mesh3D < handle
    properties
        N   % number of DOFs
        p   % node positions N by 3
        pSurface %node positions of the mesh skin
        SurfaceWeights %weights of the tets on the mesh skin, set as [] to deactivate the feature
        Surface2Tet % array of size #surface vertices storing corresponding tet for coloring the vertices | must be either used for all meshes or none
        SurfaceFaces
        surfaceVertexBaseColor %color of individual surface vertices (pre-rigidification)
        prevp
        p0  % initial node positions
        v   % node velocities
        v0  % initial node velocities
        f   % force accumulator for nodes
        triNormals %used to cache current triangle normals
        
        t      % elements vertice ids, for triangles the last entry is 0
        tetIDs % list of elements that are tets
        triIDs % list of elements that are tris
        faces  % cell array of boundaries (for collision detection)
        % Every mesh has at least one set of boundary edges
        % Each edge is a line, first column is the index of the first point
        % second column is the index of the sedgeElementsecond point
        % NOTE: order is important!!!  Treating edges as directed, the mesh
        % is to the left andt the outside is to the right, e.g., CCW order
        % for convex polygons.
        edges % vertex ids of each edge as rows edges*2 [vertex1, vertex2]. No duplicate.
        edgesRest
        edgeFaces
        facesSets % which sets of boundaries can be checked with each other for collision
        facesSetsOffset
        facesTriList
        objectDOFs %cell array of lists that containts the dofs of each individual objects
        objectDOFsPinned %same, but with pinned DOFs removed
        objectElements %cell array of per object element ids with respect to the total number of elements
        mergeMeshVerticeCount
        edgeRestDihedral %for tri3D meshes: rest angle to compute bending energy
        restEdgeAverageHeight %for tri3D meshes: rest average height to compute bending energy
        edgeNormalConsistentVertices = zeros(0,2); %for tri3D meshes : those are the vertices on the bendable edges. It is needed for normal consistency
        edgeOppositeVertices = zeros(0,2); %for tri3D meshes: opposite vertex of edge to compute bending energy
        edgeElements %for tri3D meshes:face id of the faces adjacent to the edge second entry is -1 if it is a boundary edge
        bendingEdgeRestLength
        grinspunEnergyEdgeRestLength
        bendingEdges  = zeros(0,2); % edges without the border edges
        bendingEdgesElements  = zeros(0,2); %elements of the bendable edges
        bendingEdgeAlpha1 % quick per edge alpha 1 access from average of adjacent elements alpha1
        bendingEdgeInnerAngleOrder %angle of the inner angles for bending edge when using the internalangles gptoolbox function
        bendingStiffness;
        isBendingEdge = zeros(0,2); % allows us to locate where the bending edges are
        tri3DBendingStiffness; % = 0.02;%bending stiffness for Tri3D meshes 0.5 cloth 1 paper 2 '''thin metal
        kdShellScalingFactor % bend stiffness * bendingEdgeRestLength / restEdgeAverageHeight
        prevDihedralAngle %previous dihedral angle
        prevCurvature
        restCurvature
        referenceSpaceNormals % tri3D per triangle normals
        dphidx %can be used to access individual element parts of the B matrix
        Dm     %distance matrix at rest used for strain limiting
        DmInv2D
        elementType
        vertexType

        cacheMD5 % MD5 solely use to detect change in merged meshes scene creation
        cacheMergeMaterials %solely use to detect change in merged meshes scene creation

        Graph                   % Tet adjacency graph (across edges)
        AdjagencyMatrix         % Adjacency matrix used to build the graph (not used otherwise?)
        VertexAdjacencyMatrix   % used for the quicksolve filter
        TetsPerParticle         % cellarray of the neighboring tets of each individual vertex
        
        pinned       % flags for pinned vertex indices
        pinnedInds   % indices of pinned vertices
        pinnedDOFs   % pinned indices in the 3 * N state vector
        unpinnedDOFs % unpinned indices in the 3 * N state vector
        pinnedTets   % pinned tets (tets with at least a pinned particle)
        isTetPinned  % flags for pinned tets 
        stablePinnedElements %tris with two pins or tets with 3 pins
        animationDOFs = [];
        animationInds = [];
        
        activeDOFs % REDUCED SPACE active DOFs... be really careful when using this with rigidification, go with activeDofsCorrespondingID instead... Should give the same behavior for fully elastic meshes though
        ActiveDofsCorrespondingID % list of the activeDofs positions for easier indexing
        ActiveBRows % indices of rows in B which are elastic (active)
        TriFRows %indices of triangle rows for F
        
        
        mass    % diag of mass matrix (likely more convenient?)
        M       % sparse mass matrix, 3 * N by 3 * N
        Mii
        Minv
        alpha0mass     % diag of lumped Rayleigh alpha0 multiplied with lumped mass
        Md             % sparse matrix version of alpha0mass     
        Mdii
        
        B              % kinematic relationship between vertices and element deformation gradients F
        BnoNormal      % B matrix without the normal component
        Bii
        Bt
        
        lagrangeMults  % lambdas for the compliant feedback constraint version of elasticity

        %Cf   % symbolic d2psid2F function - deprectated use mCSTVK
        restFrame   % per element material frame
        area        % per element volume without extra values like thickness
        elV         % quick per element volume access for C computation
        elMu        % quick per element Lame access for C computation
        elLambda    % quick per element lambda access for C computation
        elAlpha1    % quick per element alpha 1 access
        elBendingAlpha1
        elThickness % quick per element thickness
        elkl        % quick per element grinspun shearing stiffness
        elka        % quick per element grinspun area stiffness
        elUpper
        elLower
        vertexContactThickness

        bigAlpha1   % sparse matrix for combining with full C to build Kd
        
        RenderPatch        % plot of the elements
        FaceLines      % plot of the boundaries
        renderOffset = [0,0,0];
        
        materials           % material per tet
        materialIndex       % tet index to material
        faceMaterialIndex
        
        valence %valence of the vertices
        lineHandles %debugging line handles kept for a more efficient plotting of lines
        elementTypeEnum %an enum to store the per element type shell or tet for now
        usesTri2Tet %boolean that makes the triangles work as tets with fake 4th vertex
    end
    
    methods
        function obj = Mesh3D(p, t, attributes, materials, boundaryFaces, boundaryTIndices, SurfaceWeights, SurfaceFaces, pSurface, Surface2Tet, usesTri2Tet)
            % MESH Constructs a mesh from given points and tets.
            %   Prepares for simulation with provided materials.
            %
            %   Note that if you use meshes saved in a .mat file, you will
            %   need to regenerate the mesh if new fields get added
            %   
            if nargin == 2 || (nargin >= 3 && size(attributes,2) <= 0)
                attributes = ones(size(t,1),1);
            end
            
            if nargin <= 3 || (nargin >= 4 && size(materials,2) <= 0)
                defaultMaterial = TriangleMaterial();
            else
                defaultMaterial = materials;
            end
            
            if nargin >= 2 && size(attributes,2) > size(materials,2)
                ids = setdiff(1:size(attributes,2),1:size(materials,2));
                defaultMaterial(ids) = TriangleMaterial();
            end
            
            if isa(p, 'Mesh3D')
               % copy constructor
                fns = properties(p);
                for i = 1:numel(fns)
                    obj.(fns{i}) = p.(fns{i});
                end
                return;
            end
            
            if nargin < 5
                [boundaryFaces,obj.facesSets, boundaryTIndices] = makeBoundaries3D(p, t);
            end
            
            if nargin < 7
                SurfaceWeights = [];
                SurfaceFaces = [];
                pSurface = [];
                Surface2Tet = [];
            end

            if nargin < 11
                usesTri2Tet = false;
            end

            %creating enums
            obj.elementTypeEnum = struct('Tetrahedron',1,'Shell',2);
             if size(t,2) == 3 
                obj.elementType = t(:,1)*0 + obj.elementTypeEnum.Shell;
                obj.vertexType = zeros(obj.N,1) + obj.elementTypeEnum.Shell;
            else
                obj.elementType = t(:,1)*0 + obj.elementTypeEnum.Tetrahedron;
                obj.vertexType = zeros(obj.N,1) + obj.elementTypeEnum.Tetrahedron;
            end

            %setting mesh properties
            N = size(p, 1);
            obj.N = N;                 % number of nodes
            if usesTri2Tet
%                 BC = mexBarycenters(p,t);
%                 vert4 = BC+s*triNormals;
%                 vert4Inds = [1:size(t,1)]+N;
%                 t = [t,vert4Inds];
%                 p = [p;vert4];
            end

            obj.objectDOFs = {1:N*3};
            obj.objectDOFsPinned = obj.objectDOFs;
            obj.p = reshape(p', N * 3, 1); % position state
            obj.prevp = obj.p;
            obj.p0 = obj.p;           % initial positions
            obj.v = zeros(N * 3, 1);      % velocity state
            obj.v0 = zeros(N * 3, 1);     % initial velocities
            obj.f = zeros(N * 3, 1);      % force accumulator
            obj.mergeMeshVerticeCount = [obj.N];

            obj.SurfaceWeights = SurfaceWeights;
            obj.SurfaceFaces = SurfaceFaces;
            obj.pSurface = pSurface;
            obj.Surface2Tet = Surface2Tet;
            obj.t = int32(t);                 % tetrahedrons
            obj.objectElements = {1:size(obj.t,1)};
            obj.VertexAdjacencyMatrix = adjacency_matrix(obj.t);

            p0x1 = p(obj.t(:,1),:);
            p0x2 = p(obj.t(:,2),:);
            p0x3 = p(obj.t(:,3),:);
            normebar1 = vecnorm(p0x1-p0x2,2,2);
            normebar2 = vecnorm(p0x1-p0x3,2,2);
            normebar3 = vecnorm(p0x2-p0x3,2,2);
            obj.grinspunEnergyEdgeRestLength = [normebar1,normebar2,normebar3];

            logicalTris = obj.t(:,3) == 0;
            obj.triIDs = find(logicalTris);
            obj.tetIDs = find(~logicalTris);

            V = reshape(obj.p0, 3, size(obj.p0, 1) / 3)';
            [obj.restFrame, obj.area] = obj.makeElements3D(V, t);
            if size(obj.t,2) == 3
                obj.referenceSpaceNormals = normals(V,t,"Stable",true); 
                %gptoolbox does not normalize normals so we do it here
                obj.referenceSpaceNormals = normr(obj.referenceSpaceNormals);
            else
                obj.referenceSpaceNormals = repmat([0,0,0],size(obj.t,1),1);
            end
            
            obj.facesTriList = boundaryTIndices;
            obj.faces = boundaryFaces;
            obj.facesSets = {boundaryFaces};
            obj.facesSetsOffset = 0;
            
            obj.edges = edges(obj.faces);
            obj.edgesRest = vecnorm(p(obj.edges(:,1),:)-p(obj.edges(:,2),:),2,2);
            

            %finds the bending edges
            if size(obj.t,2) == 3
                obj.edgeElements = edge_triangle_adjacency(double(obj.t),obj.edges);
                obj.isBendingEdge = obj.edgeElements(:,2) ~= -1;
                obj.bendingEdges = obj.edges(obj.isBendingEdge,:);
                obj.bendingEdgesElements = obj.edgeElements(obj.isBendingEdge,:);
            end
            
            %computing valence
            obj.valence = zeros(obj.N,1);
            for i = 1:size(obj.t,1)
                for j = 1:size(obj.t,2)
                    vertex = obj.t(i,j);
                    obj.valence(vertex) = obj.valence(vertex) + 1;
                end
            end
            
            [obj.AdjagencyMatrix,obj.TetsPerParticle,obj.Graph] = elementAdjacencyMatrix(obj.t,N);
            
            obj.pinned = zeros(N, 1);   % flags pinned indices
            obj.pinnedInds = [];        % list of pinned node IDs
            obj.pinnedDOFs = [];
            obj.unpinnedDOFs  = 1:obj.N*3;
            obj.pinnedTets = [];
            obj.isTetPinned = zeros(size(obj.t,1),1);
            obj.activeDOFs = 1:3 * N;
            obj.ActiveDofsCorrespondingID = true(obj.N*3,1);
            obj.ActiveBRows = 1:size(t, 1) * 9;
            triNumIDs = find(obj.elementType == obj.elementTypeEnum.Shell);
            obj.TriFRows = reshape([triNumIDs*9-8,triNumIDs*9-7,triNumIDs*9-6,triNumIDs*9-5,triNumIDs*9-4,triNumIDs*9-3,triNumIDs*9-2,triNumIDs*9-1,triNumIDs*9]',1,[])';

            if size(obj.t,2) == 3 %assumes each objects has only one type of element (no tri + tet without a merge)
                [obj.Dm, obj.DmInv2D, obj.B, obj.dphidx] = computeBTri3D(V, obj.t, obj.referenceSpaceNormals);
            else
                [obj.B, obj.dphidx, obj.Dm] = computeB3D(V, obj.t);
            end

            obj.Bii = obj.B(:,obj.unpinnedDOFs);
            obj.Bt = sparse(obj.B');
            obj.BnoNormal = obj.B;

            obj.lagrangeMults = zeros(size(obj.B, 1), 1);
            
            obj.updateMaterials(attributes,defaultMaterial);
            
            if size(obj.t,2) == 3
                [obj.edgeNormalConsistentVertices, obj.edgeOppositeVertices, obj.bendingStiffness, obj.bendingEdgeInnerAngleOrder] = obj.setupMeshBending();
                [obj.edgeRestDihedral,obj.restEdgeAverageHeight, obj.bendingEdgeRestLength, obj.kdShellScalingFactor, obj.restCurvature, obj.bendingEdgeAlpha1] = obj.computeMeshBending();
                obj.prevCurvature = obj.restCurvature;
                obj.prevDihedralAngle = obj.edgeRestDihedral;
                obj.t = [obj.t,zeros(size(obj.t,1),1)];
            end

            alpha1s = reshape( repmat( obj.elAlpha1, 9, 1 ), [], 1 );
            obj.bigAlpha1 = sparse( 1:numel(alpha1s), 1:numel(alpha1s), alpha1s );
            
            obj.RenderPatch = 0;
            obj.FaceLines = {};
        end
        
        function mergeMesh(obj, mesh) 
            % MergeMesh adds a given mesh to this mesh
            % this is not necessarily efficient in that it will recompute
            % things that were computed before, but it all precomputation
            % so perhaps we don't need to worry so much.
            
            N1 = obj.N;
            N2 = mesh.N;
            numTets1 = size(obj.t,1);
            
            obj.N = N1 + N2;  
            obj.objectDOFs = {obj.objectDOFs{:}, mesh.objectDOFs{:} + (N1*3)};
            obj.objectDOFsPinned = obj.objectDOFs;
            obj.objectElements = {obj.objectElements{:}, mesh.objectElements{:} + size(obj.t,1)};
            obj.p = [ obj.p; mesh.p ];
            obj.prevp = obj.p;
            obj.p0 = [ obj.p0; mesh.p0 ];
            obj.v = [ obj.v; mesh.v ];
            obj.v0 = [ obj.v0; mesh.v0 ];
            obj.f = [ obj.f; mesh.f ];
            obj.mergeMeshVerticeCount = [obj.mergeMeshVerticeCount, mesh.mergeMeshVerticeCount];
            % triangle indices of the added mesh are offset by the number
            % of particles in this mesh
            obj.t = [ obj.t; mesh.t + N1 ];  % tets
            obj.Surface2Tet = [obj.Surface2Tet, mesh.Surface2Tet + numTets1];
            obj.elementType = [obj.elementType;mesh.elementType];
            obj.vertexType = [obj.vertexType;mesh.vertexType];
            obj.VertexAdjacencyMatrix = blkdiag(obj.VertexAdjacencyMatrix,mesh.VertexAdjacencyMatrix);

            obj.SurfaceWeights = sparse(blkdiag(obj.SurfaceWeights,mesh.SurfaceWeights));
            nFaces = numel(obj.pSurface)/3;
            obj.SurfaceFaces = [obj.SurfaceFaces;nFaces+mesh.SurfaceFaces];
            obj.pSurface = [obj.pSurface;mesh.pSurface];
            
            obj.valence = [obj.valence; mesh.valence];
            % TODO: don't really need to remake them, see fix below:
            % easiest because the elements stuipdly contain t... should 
            % refactor so that the mesh simply has Bm and area and 
            % make elements is simply a private method or part of the
            % constructor.

            % silly to recompute, but whatever... 
            V = reshape(obj.p0, 3, size(obj.p0, 1) / 3)';
            obj.area = [obj.area;mesh.area];
            obj.restFrame = cat(3,obj.restFrame,mesh.restFrame);
   
            numBoundarySets1 = numel(obj.facesSets);
            for j = 1:numel(mesh.facesSets)
                obj.facesSets{numBoundarySets1+j} = mesh.facesSets{j};
                obj.facesSetsOffset(numBoundarySets1+j) = obj.mergeMeshVerticeCount(numBoundarySets1+j-1) + obj.facesSetsOffset(numBoundarySets1+j-1);
            end
            obj.faces = [obj.faces;mesh.faces + N1]; 
            obj.facesTriList = [obj.facesTriList;mesh.facesTriList + numTets1]; 

            % concat the tet adjacency matrix and graph
            obj.AdjagencyMatrix = blkdiag(obj.AdjagencyMatrix,mesh.AdjagencyMatrix);
            obj.TetsPerParticle = [obj.TetsPerParticle;mesh.TetsPerParticle];
                       
            obj.activeDOFs = 1:3 * obj.N;
            obj.ActiveDofsCorrespondingID = [obj.ActiveDofsCorrespondingID;mesh.ActiveDofsCorrespondingID];
            obj.ActiveBRows = 1:size(obj.t, 1) * 9;
            triNumIDs = find(obj.elementType == obj.elementTypeEnum.Shell);
            obj.TriFRows = reshape([triNumIDs*9-8,triNumIDs*9-7,triNumIDs*9-6,triNumIDs*9-5,triNumIDs*9-4,triNumIDs*9-3,triNumIDs*9-2,triNumIDs*9-1,triNumIDs*9]',1,[])';
            
            obj.B = blkdiag(obj.B,mesh.B);
            V = reshape(obj.p, 3, size(obj.p, 1) / 3)';
            obj.Bt = obj.B';
            obj.BnoNormal = obj.B;
            obj.Bii = obj.B(:,obj.unpinnedDOFs);
            obj.dphidx = cat(3,obj.dphidx,mesh.dphidx);
            obj.DmInv2D = cat(3,obj.DmInv2D,mesh.DmInv2D);

            obj.lagrangeMults = zeros(size(obj.B, 1), 1);

            numEdges = size(obj.edges,1);
            obj.edges = [obj.edges;mesh.edges+numEdges];
            obj.edgesRest = [obj.edgesRest;mesh.edgesRest];
            
            %tri3D merge
            obj.edgeElements = [obj.edgeElements; mesh.edgeElements + nFaces];
            obj.isBendingEdge = logical([obj.isBendingEdge; mesh.isBendingEdge]);
            obj.bendingEdges = obj.edges(obj.isBendingEdge,:);
            obj.bendingEdgesElements = obj.edgeElements(obj.isBendingEdge,:);
            obj.edgeRestDihedral = [obj.edgeRestDihedral; mesh.edgeRestDihedral];
            obj.prevDihedralAngle = obj.edgeRestDihedral;
            obj.edgeNormalConsistentVertices = [obj.edgeNormalConsistentVertices; mesh.edgeNormalConsistentVertices + N1];
            obj.restEdgeAverageHeight = [obj.restEdgeAverageHeight; mesh.restEdgeAverageHeight];
            obj.edgeOppositeVertices = [obj.edgeOppositeVertices; mesh.edgeOppositeVertices + N1];
            obj.bendingStiffness = [obj.bendingStiffness; mesh.bendingStiffness];
            newBendingEdge = mesh.bendingEdgeInnerAngleOrder;
            newBendingEdge(:,1:2:end) = newBendingEdge(:,1:2:end) + numTets1;
            obj.bendingEdgeInnerAngleOrder = [obj.bendingEdgeInnerAngleOrder; newBendingEdge];
            obj.bendingEdgeRestLength = [obj.bendingEdgeRestLength; mesh.bendingEdgeRestLength];
            obj.grinspunEnergyEdgeRestLength = [obj.grinspunEnergyEdgeRestLength; mesh.grinspunEnergyEdgeRestLength];
            obj.kdShellScalingFactor = [obj.kdShellScalingFactor; mesh.kdShellScalingFactor];
            obj.restCurvature = [obj.restCurvature; mesh.restCurvature];
            obj.bendingEdgeAlpha1 = [obj.bendingEdgeAlpha1;mesh.bendingEdgeAlpha1];
            obj.prevCurvature = obj.restCurvature;
            
            obj.referenceSpaceNormals = [obj.referenceSpaceNormals;mesh.referenceSpaceNormals];
        
            numMat1 = numel(obj.materials);
            obj.updateMaterials([ obj.materialIndex; mesh.materialIndex + numMat1 ], [ obj.materials, mesh.materials ])
            
            alpha1s = reshape( repmat( obj.elAlpha1, 9, 1 ), [], 1 );
            obj.bigAlpha1 = sparse( 1:numel(alpha1s), 1:numel(alpha1s), alpha1s );

            obj.RenderPatch = 0;
            obj.FaceLines = {};
            obj.pin(mesh.pinnedInds + N1);
        end
        
        function clone = clone(obj)
            clone = Mesh3D(obj);
        end
        
        function updateMaterials( obj, newAttributes, newMaterials )
            % UPDATEMATERIALS Allows for new materials to be assigned to a
            % mesh after loading.
            obj.materials = newMaterials;
            obj.materialIndex = newAttributes;
            obj.faceMaterialIndex = obj.materialIndex(obj.facesTriList); 
            obj.elMu = [ obj.materials(obj.materialIndex(:)).mu ];
            obj.elkl = [ obj.materials(obj.materialIndex(:)).grinspunShearingStiffness ];
            obj.elka = [ obj.materials(obj.materialIndex(:)).grinspunAreaStiffness ];
            obj.elLambda = [ obj.materials(obj.materialIndex(:)).lambda ];
            obj.elAlpha1 = [ obj.materials(obj.materialIndex(:)).alpha1 ];
            obj.elBendingAlpha1 = [ obj.materials(obj.materialIndex(:)).bendingAlpha1 ];
            obj.elThickness = [ obj.materials(obj.materialIndex(:)).thickness ]';
            obj.elUpper = [ obj.materials(obj.materialIndex(:)).strainUpperBound ]';
            obj.elLower = [ obj.materials(obj.materialIndex(:)).strainLowerBound ]';
            isNotShell = obj.elementType ~= obj.elementTypeEnum.Shell;
            obj.elThickness(isNotShell) = 1;

            obj.vertexContactThickness = zeros(obj.N,1);
            obj.vertexContactThickness(obj.t(:,1)) = max(obj.vertexContactThickness(obj.t(:,1))  , obj.elThickness);
            obj.vertexContactThickness(obj.t(:,2)) = max(obj.vertexContactThickness(obj.t(:,2))  , obj.elThickness);
            obj.vertexContactThickness(obj.t(:,3)) = max(obj.vertexContactThickness(obj.t(:,3))  , obj.elThickness);
            obj.vertexContactThickness = obj.vertexContactThickness./2;%half thickness for contacts
            isNotShellVertex = obj.vertexType ~= obj.elementTypeEnum.Shell;
            obj.vertexContactThickness(isNotShellVertex) = 0;

            obj.elV = [obj.area]';
            obj.elV = obj.elThickness' .* obj.elV;
            
            updateMass(obj);
            obj.M = sparse( 1:3*obj.N, 1:3*obj.N, obj.mass ); % sparse matrix form for convenience
            obj.Mii = obj.M(obj.unpinnedDOFs,obj.unpinnedDOFs);
            obj.Minv = sparse( 1:3*obj.N, 1:3*obj.N, obj.mass.^-1 );
            updateAlpha0(obj);
            obj.Md = sparse( 1:3*obj.N, 1:3*obj.N, obj.alpha0mass); % sparse matrix form for convenience
            obj.Mdii = obj.Md(obj.unpinnedDOFs,obj.unpinnedDOFs);
            
            obj.surfaceVertexBaseColor = [];
            if ~isempty(obj.SurfaceWeights)
                obj.surfaceVertexBaseColor = zeros(size(obj.Surface2Tet,1),3); 
                for i=1:numel(obj.materialIndex(obj.Surface2Tet))
                    tetID = obj.Surface2Tet(i);
                    matID = obj.materialIndex(tetID);
                    obj.surfaceVertexBaseColor(i,:) = obj.materials(matID).color;
                end
            end
        end
            
        function x = getPositionFormatted(obj)
            %Get x in a n x 3 format compatible with most gptoolbox
            %functions
            x = reshape(obj.p',3,obj.N)';
        end

        function x = getPosition0Formatted(obj)
            %Get x in a n x 3 format compatible with most gptoolbox
            %functions
            x = reshape(obj.p0',3,obj.N)';
        end

        function [B, gamma] = getB(obj, ~)
           B = obj.B;
           gamma = speye(obj.N * 3);
        end

        function [B, gamma] = peekB(obj, cache, bodyPositions, meshp)
           [B, gamma] = obj.getB();
        end
        
        function M = getM(mesh)
            M = mesh.M;
        end
        
        function f = getCurrentForce(mesh)
            f = mesh.f;
        end
        
        function v = getCurrentVelocity(mesh)
            v = mesh.v;
        end
        
        function linearMomentum = getLinearMomentum( mesh )
            linearMomentum = sum( reshape(mesh.mass.*mesh.v, 3, mesh.N), 2 );
        end
        
        function labelTriangles( mesh )
             pt = reshape( mesh.p, 3, mesh.N );
             for i=1:size(mesh.t,1)
                 pos = sum( pt( :, mesh.t(i,:) ), 3 ) / 4;
                 text( pos(1), pos(2), pos(3), ""+i );
             end        
        end

        function labelVertices( mesh )
             pt = reshape( mesh.p, 3, mesh.N );
             for i=1:mesh.N
                 pos = pt(:,i);
                 text( pos(1), pos(2), pos(3), ""+i );
             end        
        end

        function updateMass(mesh)
            % UPDATEMASS is used to contsruct the mass damping and is
            % typically only called on construction of the mesh
            mesh.mass = computeMass(mesh);
        end
        
        function updateAlpha0( mesh )
            % UPDATEALPHA0 is used to contsruct the mass damping and is
            % typically only called on construction of the mesh
            
            alpha0tet = [mesh.materials(mesh.materialIndex).alpha0];
            %values are declared per tet, but we want them per node
            alpha0vector = zeros( mesh.N*3, 1 );
            for i = 1:size(mesh.t,1)
                dim = 4;
                if mesh.elementType(i) == mesh.elementTypeEnum.Shell
                    dim = 3;
                end
                for node = 1:dim
                    id = mesh.t(i,node);
                    value = alpha0vector(id*3) + 1/4 * alpha0tet(i);
                    alpha0vector(id*3) = value;
                    alpha0vector(id*3-1) = value;
                    alpha0vector(id*3-2) = value;
                end
            end
            mesh.alpha0mass = alpha0vector .* mesh.mass ;
        end
        
        function [edgeNormalConsistentVertices, edgeOppositeVertices,bendingStiffness, bendingEdgeInnerAngleOrder] = setupMeshBending(mesh3D)
            %Discrete Shells Eitan Grinspun Eurographics/SIGGRAPH Symposium on Computer Animation
            
            edgeNormalConsistentVertices = zeros(size(mesh3D.bendingEdges,1),2);
            edgeOppositeVertices = zeros(size(mesh3D.bendingEdges,1),2);
            bendingStiffness = zeros(size(mesh3D.bendingEdges,1),1);
            bendingEdgeInnerAngleOrder = zeros(size(mesh3D.bendingEdges,1),8);
            counter = 1;

            for i = 1:size(mesh3D.bendingEdges,1)
                vertex1 = mesh3D.bendingEdges(i,1);
                vertex2 = mesh3D.bendingEdges(i,2);
                trisOfEdge = mesh3D.bendingEdgesElements(i,:);
                
                v1 = mesh3D.p(vertex1*3-2:vertex1*3);
                v2 = mesh3D.p(vertex2*3-2:vertex2*3);

                tri1 = mesh3D.t(trisOfEdge(1),:);
                vertexOpposite1 = tri1 ~= vertex1 & tri1 ~= vertex2;
                vertexOpposite1 = tri1(vertexOpposite1);

                tri2 = mesh3D.t(trisOfEdge(2),:);
                vertexOpposite2 = tri2 ~= vertex1 & tri2 ~= vertex2;
                vertexOpposite2 = tri2(vertexOpposite2);

                o1 = mesh3D.p(vertexOpposite1*3-2:vertexOpposite1*3);
                o2 = mesh3D.p(vertexOpposite2*3-2:vertexOpposite2*3);

                if ~mesh3D.isEdgeNormalsConsistent(o1,v1,v2,o2, mesh3D.referenceSpaceNormals(trisOfEdge(1),:),mesh3D.referenceSpaceNormals(trisOfEdge(2),:))
                    %flip the edge, because it is inconsistent with the
                    %mesh normals
                    vertexTmp = vertex1;
                    vertex1 = vertex2;
                    vertex2 = vertexTmp;
                end
                edgeNormalConsistentVertices(counter,1) = vertex1;
                edgeNormalConsistentVertices(counter,2) = vertex2;
                edgeOppositeVertices(counter,1) = vertexOpposite1;
                edgeOppositeVertices(counter,2) = vertexOpposite2;

                %compute the bending stiffness as an average of its
                %triangle materials
                mat1 = mesh3D.materials(mesh3D.materialIndex(trisOfEdge(1)));
                mat2 = mesh3D.materials(mesh3D.materialIndex(trisOfEdge(2)));
                [E1,nu1] = toEmu(mat1.lambda,mat1.mu);
                [E2,nu2] = toEmu(mat2.lambda,mat2.mu);
                E = (E1 + E2)/2;
                nu = (nu1 + nu2)/2;
                thickness = (mat1.thickness + mat2.thickness)/2;
                bendingStiffness(i) = materialToBendingStiffness(E,nu,thickness);

                for j = 1:3
                    if tri1(j) == vertex1
                        bendingEdgeInnerAngleOrder(i,1) = trisOfEdge(1);
                        bendingEdgeInnerAngleOrder(i,2) = j;
                    end
                    if tri2(j) == vertex1
                        bendingEdgeInnerAngleOrder(i,3) = trisOfEdge(2);
                        bendingEdgeInnerAngleOrder(i,4) = j;
                    end
    
                    if tri1(j) == vertex2
                        bendingEdgeInnerAngleOrder(i,5) = trisOfEdge(1);
                        bendingEdgeInnerAngleOrder(i,6) = j;
                    end
                    if tri2(j) == vertex2
                        bendingEdgeInnerAngleOrder(i,7) = trisOfEdge(2);
                        bendingEdgeInnerAngleOrder(i,8) = j;
                    end
                end

                counter = counter + 1;
            end
        end

        function [edgeRestDihedral,restDualGraphEdgeLength, bendingEdgeRestLength, kdShellScalingFactor, curvature, bendingEdgeAlpha1]=computeMeshBending(mesh3D)
            [curvature, edgeRestDihedral, BC, restDualGraphEdgeLength, bendingEdgeRestLength]= computeCurvature(mesh3D,mesh3D.p);
            %https://la.disneyresearch.com/wp-content/uploads/Discrete-Bending-Forces-and-Their-Jacobians-Paper.pdf
            %page 5 section 5
            restAreaSum = (mesh3D.area(mesh3D.bendingEdgesElements(:,1)) + mesh3D.area(mesh3D.bendingEdgesElements(:,2)));
            restScaling = (3.*bendingEdgeRestLength.*bendingEdgeRestLength) ./ restAreaSum;
            kdShellScalingFactor = mesh3D.bendingStiffness .* restScaling;
            bendingEdgeAlpha1 = sum(mesh3D.elBendingAlpha1(mesh3D.bendingEdgesElements(:,:)),2)./2;
        end
        
        %% public function prototypes
        prepare( obj, infMassForPinned )
        pin( obj, pinnedInds )
        updateParticles( obj, h, deltav )
        applyAcceleration( obj, acc )
        resetForce( obj )
        setRigidTransform( obj, degrees, translation, scale )
        setRigidMotion( obj, omega, radPerSec, velocity  )
        setRigidRotationFromMatrix(obj, rotationMatrix, forceCenterOfMass)
        computeActiveDOFs(obj) 
        f = peekCurrentForce(mesh, h, AngularVelocity, RigidBodyPosition, Inertia, meshf, meshp);
    end
    methods(Static)
        [W, T, SV, V, V2Htet] = voxelizer(V,FACES,side);
        [W, V2Htet] = meshSkinning(Vsurface,TetMesh3D);
        [restFrame, area] = makeElements3D(V, T);
        

        function isConsistent = isEdgeNormalsConsistent(x0, x1,x2,x3,normal1,normal2)
            e0 = (x2-x1);
            e2 = (x0-x1);
            eTilde2 = (x3-x1);

            [edgeNormal1, edgeNormal2] = computePerEdgeFaceNormals(e0',e2',eTilde2');
            check1 = norm(normal1 - edgeNormal1) < 0.000000001;
            check2 = norm(normal2 - edgeNormal2) < 0.000000001;
            isConsistent = check1 && check2;
        end

        function x = formatPositions(p)
            %Get x in a n x 3 format compatible with most gptoolbox
            %functions
            x = reshape(p',3,[])';
        end
    end
end

