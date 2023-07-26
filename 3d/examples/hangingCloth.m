clear;
cla;
clear meshes;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

rho = 2;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.01;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness)];

baseMesh = generateCloth([], tMaterial, 0.2);
settings.recomputeCacheAinv = true;

% pinning tris
zPos = baseMesh.p(3:3:end);
max_I = find(zPos >= max(zPos)-0.01);
baseMesh.pin(sort(max_I));

baseMesh.setRigidTransform([-30,0,0],[0,0,0]);
meshes = AdaptiveMesh3D(baseMesh);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 0;
rigidificator.ElastificationThreshold = 0; 
rigidificator.setBendingThresholdsFromPlanar();

settings.RigidificationEnabled = false; % There are no joint constraints so this would stop the motion of the rigid bodies where only the top is pinned. TODO: Need to handle those

integrator = LDLBackwardEuler3D();

integrator.Gravity = -9.8;
energyModel = StVenantKirchoff3DEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.SceneName = 'hangingCloth';
settings.addShellNormalDeformation = 1;
settings.addBendingEnergy = 1; %Note: the bending energy is very slow to
% compute with matlab. Using a mex version would speed it up.
settings.campos=[5,5,0];

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);