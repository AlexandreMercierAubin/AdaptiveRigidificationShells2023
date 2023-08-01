clear;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear meshes;
rho = 3;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.005;   % Rayleigh factor on M
alpha1 = 0.1;  % Rayleigh factor on K
thickness = 0.015;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
ka = 0.1;
kl = 0.15;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8],strainUpperBound,strainLowerBound,thickness,kl,ka,0.1)];

resetMesh = true;
baseMesh = shellOBJLoader('bunnyLow',[],tMaterial,[1,1,1],resetMesh,settings);
% generateCloth([], tMaterial, 0.1);
baseMesh.setRigidTransform([90,0,0],[0,0,0.3],true);

% pinning tris

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDotClothRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 1e-3; 
rigidificator.setBendingThresholdsFromPlanar();
integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0,0);
integrator.useQuicksolveContactFilter = 4;

energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,-0.2], 0.7);
contactFinder = {planeContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 15/h;
% settings.PlotEDotHist = 1;
settings.SceneName = 'clothBunnyGround';
settings.WriteOBJs = true;
settings.OBJDir = './objs/clothBunnyGround/';
settings.campos=[5,5,.5];
settings.StrainLimitingEnabled = true;
settings.addBendingEnergy = 1;
settings.PCGiterations = 1;
settings.useGrinspunPlanarEnergy = true;

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);