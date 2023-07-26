cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear mesh3Da;

rho = 10;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
k = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, k);
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
strainUpperBound = 1.3;
strainLowerBound = 0.7;
thickness = 0.1;
ka = 5;
kl = 50; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka)];


V = [1,-1,0.5;
     1,0,0;
     1,0,1]*0.5;
T = [1,2,3];
J = [1]';

mesh3D = Mesh3D(V,T,[],tMaterial,T,J,[],[],[],[]);

mesh3D.setRigidTransform([90,0,-10],[-1,0,1],true);%the tet is just a hack so the scene renders
% baseMesh.pin([1,3,4]);

mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-3;
rigidificator.ElastificationThreshold = 1e-2; 
rigidificator.setBendingThresholdsFromPlanar();
% integrator = FullNewton3D();
integrator = LDLBackwardEuler3D();
% integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.001,0.00001);
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.0,0.0);
% energyModel = NeoHookean3DEnergy();
energyModel = StVenantKirchoff3DEnergy();

mannequinPosition = [0,0,0];
plotRatio= 0.99;
frictionCoeff=0.8;
animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/platformSine.mot','./3d/data/AnimationData/platform/',[90,0,0],[0.5,0.5,0.5],mannequinPosition,frictionCoeff,h);
animContactFinder.plotScale = plotRatio;
% animContactFinder.stopFrame = 1;

animContactFinder.interpolateOnFrames = 5;
% animContactFinder.skipdetection = true;
meshMeshContactFinder = MeshSCD3D(0.0,false, []);
contactFinder = {animContactFinder,meshMeshContactFinder};% 

settings.MakeVideo = 1;
settings.SceneName = 'animatedCollider';
% settings.FramesToRecord = 300;
% settings.PCGiterations = 100;
% settings.quicksolveSimulation = true;
settings.campos=[5,5,1];
settings.camtarget = [0,0,1];
% settings.RigidificationEnabled = false;
settings.addShellNormalDeformation = false;
settings.recomputeCacheAinv = true;
settings.useGrinspunPlanarEnergy = true;
% integrator.separateQuicksolveGravity = false;
% settings.quicksolveSimulation = true;

td = simulate3D({mesh3Da},h,contactFinder, integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"singleTet",["default","adaptive"]);