cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 3;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.1;
strainLowerBound = 0.9;
ka = 5;
kl = 2; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka)];

resetMesh = false;
scaling = 0.98;

nSides = 20;
length = 6;
initialRadius = 0.65;
sineAlpha = 0.02;
sineFreq = 0.8;
phase = (-pi/2)/sineFreq;
baseMesh = generateShellTube([],tMaterial,nSides,length,initialRadius,sineAlpha,sineFreq, phase, [1,1,1]);

baseMesh.setRigidTransform([0,-90,0],[0,0,0],false);
baseMesh.setRigidTransform([0,0,0],[2.7,0,0],true);
% baseMesh.pin([1,3,4]);

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 5e-3;
rigidificator.ElastificationThreshold = 5e-2; 
rigidificator.setBendingThresholdsFromPlanar();
% integrator = FullNewton3D();
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0, 0);
integrator.useFullAinv = true;
energyModel = NeoHookean3DEnergy();
% energyModel = StVenantKirchoff3DEnergy();

frictionCoeff = 0.1;
animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/armFlex.mot','./3d/data/AnimationData/arm/',[0,0,-90],[1.1,1,1.1],[0,0,0],frictionCoeff,h);
animContactFinder.stopFrame = 1;%stop at second 6 of the motion
animContactFinder.plotScale = 0.95;
animContactFinder.interpolateOnFrames = 0.05/h;
% animContactFinder.skipdetection = true;
contactFinder = {animContactFinder};% 

settings.MakeVideo = 1;
settings.SceneName = 'clothArmStatic';
% settings.FramesToRecord = 300000;
settings.PCGiterations = 1;
% settings.quicksolveSimulation = true;
settings.campos=[5,10,0]*2.5;
settings.camtarget = [0,0,0];
% settings.RigidificationEnabled = false;
settings.recomputeCacheAinv = true;
% settings.StrainLimitingEnabled = true;
settings.renderer = 'opengl';
% settings.PlotSkip = 10;
% settings.WriteOBJs = true;
settings.OBJDir = './objs/clothArm/';
settings.PGSiterations = 25;
% settings.PlotSkip = plotSkip60FPS(h);
settings.addBendingEnergy= true;
% settings.addShellNormalDeformation = true;
settings.useGrinspunPlanarEnergy = true;

% integrator.separateQuicksolveGravity = false;
% integrator.PGSquicksolve = true;
% settings.PCGiterations = 200;
% settings.quicksolveSimulation = true;
% integrator.useFullContactAinv = true;

td = simulate3D({meshes},h,contactFinder, integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"singleTet",["default","adaptive"]);