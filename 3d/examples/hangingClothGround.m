clear;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear meshes;
rho = 10;
nu = 0.33;      % Poisson ratio: close to 0.5 for rubber
E = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
thickness = 0.01;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness)];

baseMesh = generateCloth([], tMaterial, 0.1);

% pinning tris
zPos = baseMesh.p(3:3:end);
max_I = find(zPos >= max(zPos)-0.01);
baseMesh.pin(sort(max_I));

baseMesh.setRigidTransform([-35,0,0],[0,0,0.4],true);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDotClothRigidificator();
rigidificator.RigidificationThreshold = 1e-5;
rigidificator.ElastificationThreshold = 2e-4; 
rigidificator.RigidificationBendThreshold = 1e-1;
rigidificator.ElastificationBendThreshold = 1e-2;
rigidificator.PreventPinnedRigidification = true; %pinned dofs won't rigidified because it could lock joints (need to implement joint constraints)
integrator = BackwardEuler3D();
integrator.Gravity = -9.8;

energyModel = CorotationalEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,-0.1], 0.7);
contactFinder = {planeContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 2000;
settings.InitialWindowPosition = [0,0,1920,1080];
settings.SceneName = 'hangingClothGround';
settings.addShellNormalDeformation = 1;
settings.campos=[5,5,.5];
settings.recomputeCacheAinv = true;
settings.addBendingEnergy = true;

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);