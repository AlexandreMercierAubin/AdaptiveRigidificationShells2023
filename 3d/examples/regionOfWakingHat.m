clear
close all

h = 0.003;

% CONFIG
settings = Simulation3DSettings();
settings = Simulation3DSettings();
settings.FirstFrameRigidification = true;
settings.useGrinspunPlanarEnergy = true;
settings.campos=[2,2,7];

% MESHES
rho = 1;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0005;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
color = [0.3,0.46,0.8];
ka = 1;
kl = 5; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka)];

%TOP
rho = 2;
ka2 = 10;
kl2 = 50; 
thickness = 0.1;
nu2 = 0.4;      % Poisson ratio: close to 0.5 for rubber
E2 = 5e4;
[ mu2, lambda2 ] = toLame( nu2, E2 );
alpha0 = 0.005;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu2, lambda2, alpha0, alpha1, [0.3,0.2,0.8], strainUpperBound, strainLowerBound, thickness, kl2, ka2)];
tMaterial= [tMaterial,tMaterial2];

resetMesh = true;
mesh3D = shellOBJLoader('hathd',[],tMaterial,[0.5,0.5,0.5],resetMesh,settings);

v = reshape(mesh3D.p,3,[]);
n1 = v( :, mesh3D.t(:,1) );
n2 = v( :, mesh3D.t(:,2) );
n3 = v( :, mesh3D.t(:,3) );
elCenter = (1/3) * (n1+n2+n3);
% The .24 here is because it doesn't quite sit on the y axis... small shift
% in the z direction.
dfromyaxis = sqrt(sum((elCenter([1,3],:)-[0;0.05]).^2, 1));
dfromyaxis2 = sqrt(sum((elCenter([1,3],:)-[0;-0.1]).^2, 1));
rthresh = 0.3;
rthresh2 = 0.34;
% find the bits away from the core, and avoid the boxes on the side of teh
% rocket or other features that are too far off the centerline
ind =  (dfromyaxis > rthresh  & abs(elCenter(2,:))<1) &(dfromyaxis2 > rthresh2 & abs(elCenter(2,:))<1); % this is inverted as they are doubles, not logical

attributes = mesh3D.materialIndex;
attributes(ind) = 2;
mesh3D.updateMaterials( attributes, tMaterial );
mesh3D.setRigidTransform([70,30,0],[0,0,0.4],true);
mesh3Da = AdaptiveMesh3D(mesh3D);

cache = QuickSolveCache3D();
energyModel = NeoHookean3DEnergy();

%preconditioner for quicksolve
%compute young's moduli
precomputeAInverseDiagonalBlocks3D( cache, mesh3Da, h, energyModel, settings);
G = mesh3Da.elMu;
lam = mesh3Da.elLambda;
youngs_moduli = (G.*(6.*lam+2.*G))./(lam+G);
Apre = mesh3Da.M + h.*h.*...
             mesh3Da.B'*...
             spdiags(reshape(repmat(youngs_moduli, 9,1), 9*numel(youngs_moduli),1), [0], 9*numel(youngs_moduli), 9*numel(youngs_moduli))*...
             mesh3Da.B;
cache.Apre = Apre;
%pinned vertices
ii = mesh3Da.unpinnedDOFs;

% prep the matrix and permutation
Aii = Apre(ii,ii);
perm = symamd(Aii);
App = Aii(perm,perm);
I = speye(size(Aii));
Perm = I(perm,:);

Lpre = ichol( App, struct('type','nofill','michol','on') ); 
cache.Lpre = Lpre;
cache.preconditioner = @(r) Perm'*(Lpre'\(Lpre\(Perm*r)));


isShellElement = mesh3Da.elementType == mesh3Da.elementTypeEnum.Shell;
mesh3Da.triNormals = zeros(size(mesh3D.t,1),3);
mesh3Da.triNormals(isShellElement,:) = mexNormals(mesh3Da.getPositionFormatted,mesh3Da.t(isShellElement,1:3));
[psi, cache.bendingForces,cache.shellBendingH, cache.shellBendingDampingH, cache.dihedralAngles] = computeBendingGradHess(mesh3Da, mesh3Da.p, mesh3Da.triNormals, cache, settings);       
isShell = mesh3Da.elementType == mesh3Da.elementTypeEnum.Shell;
[ ii, jj, vals, dpsidX , valsD, Wa, Wl] = mexGrinspunPlanar(mesh3Da.p,mesh3Da.t,mesh3Da.elkl,mesh3Da.elka,mesh3Da.area,isShell, mesh3Da.elAlpha1, mesh3Da.grinspunEnergyEdgeRestLength, sum(isShell));
cache.grinspunPlanarH = -sparse(ii,jj,vals,mesh3Da.N*3,mesh3Da.N*3);
cache.grinspunPlanarHd = -sparse(ii,jj,valsD,mesh3Da.N*3,mesh3Da.N*3);
cache.grinspunPlanarForces = -dpsidX;
allpsi = Wa + Wl;
sumWa = sum(Wa);
sumWl = sum(Wl);
psi = sumWa + sumWl;
cache.K = cache.grinspunPlanarH +cache.shellBendingH;
cache.Kd = cache.grinspunPlanarHd + cache.shellBendingDampingH;
cache.D = cache.Kd + mesh3D.Md;

cache.F = mexComputeFtri3D(mesh3Da.getPositionFormatted,mesh3Da.t,mesh3Da.triNormals,mesh3Da.dphidx(:));
be = LDLBackwardEuler3D();
be.useQuicksolveContactFilter = 5;
settings.PCGiterations = 1;
be.Gravity = -9.8;
% be.PGSquicksolve = true;
be.useFullContactAinv = true;
be.useFullAinv = true;

rigidificator = ECurvCloth3DRigidificator();
rigidificator.bendType = 2;
rigidificator.FrameCount = 1;


gcf;
[mainFig,axesList, initialCamera] = setupWindow(settings,{mesh3Da});
hold on;
axis off;

impulseVertIDs = [152,1101];
impulse = [0.05,0.05];
impulsePhi = 0.000; %positive phi so it deformes a tiny bit when gravity is factored in
phi = [impulsePhi;
  impulsePhi];
Jc = sparse([3;2;1;6;5;4],[454;455;456;3301;3302;3303],[-1;1;1;-1;1;1],6,4263);

% colorbarHandle.Limits = [-10 0];
axis([-2, 2, -4, 0]);
plotMesh3D(mesh3Da,cache,settings);
pr = mesh3Da.getPositionFormatted;
V1 = pr(impulseVertIDs(1),:);
xdata = [V1(1); V1(1)];
ydata = [V1(2); V1(2)];
zdata = [V1(3); V1(3)-impulse(1)];
line(xdata,ydata,zdata, 'Color',[100/255,100/255,255/255])

V2 = pr(impulseVertIDs(2),:);
xdata = [V2(1); V2(1)];
ydata = [V2(2); V2(2)];
zdata = [V2(3); V2(3)-impulse(2)];
line(xdata,ydata,zdata, 'Color',[100/255,100/255,255/255])

computeRegionOfWaking(mesh3Da, h, cache, be, settings, rigidificator, V1, V2, impulseVertIDs, impulse, Jc, phi);
cache.ApproximatedDeltaV = zeros(size(cache.ApproximatedDeltaV));


function computeRegionOfWaking(mesh, h, cache, be, settings, rigidificator, V1, V2, impulseVertIDs, impulse, Jc, phi)
    cache.clearWarmStartInfo();
    
    cInfo(1) = contactInfo3D( [V1(1), V1(2), V1(3)+impulse(1)], [0,0,1], [0,1,0], 0, impulseVertIDs(1), 1, [1,0,0]);
    cInfo(2) = contactInfo3D( [V2(1), V2(2), V2(3)+impulse(2)], [0,0,1], [0,1,0], 0, impulseVertIDs(2), 1, [1,0,0]);
    cache.cInfo = cInfo;
    quickSolve3D( cache, be, mesh, h, Jc, phi, settings, [], 1, cInfo); 

    faceColor = -8*ones(size(mesh.t,1),1);
    for i = 8:-1:0
        mesh.RigidificationValues = zeros(size(mesh.isTetElastic,1),rigidificator.FrameCount);
        
        makeMeshRigid(mesh);
        rigidificator.ElastificationThreshold = 10 ^ (-i);
        rigidificator.RigidificationThreshold = 1;
        rigidificator.setBendingThresholdsFromPlanar();
        rigidificator.Fprev = cache.F;
        rigidificator.checkForElastification(mesh, cache, rigidificator.FrameCount+1, h, settings);
        elasticTris = mesh.ElasticTetInds;
        faceColor(elasticTris) = -i;
    end

    mesh.RenderPatch.FaceVertexCData = faceColor;
    downwardColor = linspace(0,1,200)';
    upwardColor = linspace(1,0,200)';
    redChannel = [upwardColor];
    greenChannel = [downwardColor];
    customCMap = [redChannel,greenChannel, zeros(200,1)];
    colorbarHandle = colorbar;
    clim([-8,-5]);
    colormap(customCMap);
    set(colorbarHandle,'position',[.9 .3 .04 .5])
    pause(0.4);
end

function makeMeshRigid(meshTri3D)
    rigidIDbyVert = int32(ones(meshTri3D.N,1));
	[ com, comdot, mass, rotMass, angularMomentum, vertexDisp] = mexRigidBodyProperties3D( 1, rigidIDbyVert, meshTri3D.p, meshTri3D.v, meshTri3D.mass );

    rigidVert = 1:meshTri3D.N;
    inds3 = rigidVert * 3;
    meshTri3D.VertexDisp(rigidVert,:) = [vertexDisp(inds3-2), vertexDisp(inds3-1), vertexDisp(inds3)];
    numRigid =1;
    rigidIDbyElement = ones(size(meshTri3D.t,1),1);

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
        meshTri3D.RigidBodies(i).isPinned = false;
        meshTri3D.RigidBodies(i).isAnimated = false; %Todo: implement scripted animations for cloth
    end
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
    meshTri3D.ElasticTetInds = [];
    meshTri3D.ElasticInds = [];
    meshTri3D.ActiveBRows = [];
    meshTri3D.isTetElastic = false(size(meshTri3D.t,1),1);
    meshTri3D.computeActiveDOFs();
end 