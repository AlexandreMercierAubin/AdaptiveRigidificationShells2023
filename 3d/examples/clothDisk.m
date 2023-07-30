clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesh3Da;
rho = 3;
nu = 0.36;      % Poisson ratio: close to 0.5 for rubber
E = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.05;
strainUpperBound = 1.2;
strainLowerBound = 0.8;
ka = 0.01;
kl = 2.5; %under 2.5, it will break
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.3,0.3,0.7],strainUpperBound, strainLowerBound,thickness, kl, ka)];

resetMesh = true;
mesh3D = shellOBJLoader('disk',[],tMaterial,[0.5,0.5,0.5],resetMesh,settings);
mesh3D.setRigidTransform([60,0,0],[0,0,1]);

% pinning tris
mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 5e-4;
rigidificator.ElastificationThreshold = 5e-3; 
rigidificator.setBendingThresholdsFromPlanar();

integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.0001, 0.0001);


energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,-0.2], 0.2);
contactFinder = {planeContactFinder};

settings.MakeVideo = 1;
settings.SceneName = 'disk';
settings.campos=[5,-2,.75];
settings.addBendingEnergy = 1;
settings.useGrinspunPlanarEnergy = true;
settings.PGSiterations = 30;

integrator.useFullAinv = true;
integrator.useFullContactAinv = true;
integrator.useQuicksolveContactFilter = 7;
settings.recomputeCacheAinv = true;

simulate3D(mesh3Da,h,contactFinder, integrator, rigidificator, settings, energyModel);