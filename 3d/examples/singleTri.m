clear;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear meshes;
rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
strainUpperBound = 1.3;
strainLowerBound = 0.7;
thickness = 0.1;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[1,0,0.5], strainUpperBound, strainLowerBound,thickness)];

V = [1,-1,0.5;
     1,0,0;
     1,0,1];
T = [1,2,3];
J = [1]';

baseMesh = Mesh3D(V,T,[],tMaterial,T,J,[],[],[],[]);

meshes = AdaptiveMesh3D(baseMesh);
meshes.AlwaysRigidElements = [1:size(meshes.t,1)]';

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-2;
rigidificator.ElastificationThreshold = 2e-2; 
rigidificator.setBendingThresholdsFromPlanar();
integrator2 = LDLBackwardEuler3D();
integrator2.Gravity = -0.5;

energyModel = NeoHookean3DEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.RigidificationEnabled = false;
settings.MakeVideo = 1;
settings.SceneName = 'twoTris3D';
settings.addShellNormalDeformation = 1;
settings.addBendingEnergy = true;
settings.campos=[5,10,8];
settings.recomputeCacheAinv = true;
settings.DrawForces = true;

simulate3D(meshes,h,contactFinder, {integrator2}, rigidificator, settings, energyModel);