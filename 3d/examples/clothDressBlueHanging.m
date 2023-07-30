cla;
clear;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 4;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.01;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness)];

resetMesh = true;
baseMesh = shellOBJLoader('dressblueCoarse',[],tMaterial,[1,1,1],resetMesh,settings);
z = baseMesh.p(3:3:end);
vertIDs = find(z >= max(z)-0.05);
baseMesh.pin(vertIDs);

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = EDotClothRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 1e-3; 
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.1,0.1 );
energyModel = NeoHookean3DEnergy();

frictionCoeff = 0.3;

settings.MakeVideo = 1;
settings.SceneName = 'dressblueHanging';
settings.PGSiterations = 20;
settings.campos=[-3.5,-3.5,0.35];
settings.camtarget = [0,0,0.1];
settings.renderer = 'opengl';
settings.camLightPosition = 'headlight';
settings.useGrinspunPlanarEnergy = true;

td = simulate3D({baseMesh},h,NullContactFinder(), integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"singleTet",["default","adaptive"]);