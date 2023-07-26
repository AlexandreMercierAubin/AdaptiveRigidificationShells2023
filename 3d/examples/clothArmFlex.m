cla;
clear;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear meshes;

rho = 5;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.1;
strainLowerBound = 0.9;
ka = 40;
kl = 20; 
bendingAlpha1 = 0.05;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,bendingAlpha1)];

resetMesh = false;
scaling = 0.98;

nSides = 35;%This gives good wrinkles
% nSides = 15;%This is a better test resolution
length = 8;
initialRadius = 0.65;
sineAlpha = 0.00;
sineFreq = 0.8;
phase = (-pi/2)/sineFreq;
baseMesh = generateShellTube([],tMaterial,nSides,length,initialRadius,sineAlpha,sineFreq, phase, [1,1,1]);

baseMesh.setRigidTransform([0,-90,0],[0,0,0],false);
baseMesh.setRigidTransform([0,0,0],[3.1,0,0],true);
% baseMesh.pin([1,3,4]);

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = ECurvCloth3DRigidificator();
% rigidificator.RigidificationThreshold = 7e-3; %with curvatures
% rigidificator.ElastificationThreshold = 7e-1; 
rigidificator.RigidificationThreshold = 2e-3;% with angles
rigidificator.ElastificationThreshold = 2e-2; 
rigidificator.setBendingThresholdsFromPlanar();

% integrator = FullNewton3D();
integrator = LDLBackwardEuler3D();
% integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.001,0.001);
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0, 0.0);
% integrator.useFullAinv = true;
energyModel = NeoHookean3DEnergy();
% energyModel = StVenantKirchoff3DEnergy();

frictionCoeff = 0.5;
animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/armFlex.mot','./3d/data/AnimationData/armPretty/',[0,0,-90],[1.1,1,1.1],[0,0,0],frictionCoeff,h);
% animContactFinder.stopFrame = 6/h;%stop at second 6 of the motion
animContactFinder.plotScale = 0.95;
animContactFinder.interpolateOnFrames = 0.05/h;
% animContactFinder.skipdetection = true;
contactFinder = {animContactFinder};% 

settings.MakeVideo = 1;
settings.SceneName = 'clothArmFlex';
settings.FramesToRecord = 30/h;%steps to record...
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
settings.OBJDir = './objs/clothArmDefault/';
settings.PGSiterations = 20;
settings.PlotSkip = plotSkip60FPS(h);
settings.addBendingEnergy= true;
% settings.addShellNormalDeformation = true;
settings.useGrinspunPlanarEnergy = true;
rigidificator.bendType = 2;
% integrator.useQuicksolveContactFilter = 2; %special CG warmstart
% settings.PCGiterations = 100;
% integrator.useFullContactAinv = true;

td = simulate3D({baseMesh,meshes},h,contactFinder, integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"armFlex_",["default","adaptive"]);%todo rename singletet to this oopsi
save("armFlex"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
readTDcsvLog(["armFlex_default.csv","armFlex_adaptive.csv"],["blue","#FF6600"],-1,h);