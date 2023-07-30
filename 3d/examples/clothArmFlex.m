cla;
clear;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
clear mesh3Da;

rho = 5;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.1;
strainLowerBound = 0.9;
ka = 70;
kl = 40; 
bendingAlpha1 = 0.1;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,bendingAlpha1)];

resetMesh = false;
scaling = 0.98;

nSides = 35;%This gives good wrinkles
% nSides = 15;%This is a better test resolution
length = 8;
initialRadius = 0.60;
sineAlpha = 0.00;
sineFreq = 0.8;
phase = (-pi/2)/sineFreq;
mesh3D = generateShellTube([],tMaterial,nSides,length,initialRadius,sineAlpha,sineFreq, phase, [1,1,1]);

mesh3D.setRigidTransform([0,-90,0],[0,0,0],false);
mesh3D.setRigidTransform([0,0,0],[3.1,0,0],true);

mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-3;% with angles
rigidificator.ElastificationThreshold = 5e-3; 
rigidificator.ElastificationBendThreshold = 2e-1;
rigidificator.RigidificationBendThreshold = 2e-2;

integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0, 0.0);
integrator.useQuicksolveContactFilter = 4;
energyModel = NeoHookean3DEnergy();

frictionCoeff = 0.5;
animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/armFlex.mot','./3d/data/AnimationData/armPretty/',[0,0,-90],[1.1,1,1.1],[0,0,0],frictionCoeff,h);
animContactFinder.plotScale = 0.95;
animContactFinder.interpolateOnFrames = 0.05/h;
contactFinder = {animContactFinder};

settings.MakeVideo = 1;
settings.SceneName = 'clothArmFlex';
settings.FramesToRecord = 30/h;%steps to record...
settings.PCGiterations = 1;
settings.campos=[5,10,0]*2.5;
settings.camtarget = [0,0,0];
settings.recomputeCacheAinv = true;
settings.renderer = 'opengl';
settings.WriteOBJs = true;
settings.OBJDir = './objs/clothArmAdaptive/';
settings.PGSiterations = 20;
settings.PlotSkip = plotSkip60FPS(h);
settings.addBendingEnergy= true;
settings.useGrinspunPlanarEnergy = true;
rigidificator.bendType = 2;

td = simulate3D({mesh3D,mesh3Da},h,contactFinder, integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"armFlex_",["default","adaptive"]);
save("armFlex"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
readTDcsvLog(["armFlex_default.csv","armFlex_adaptive.csv"],["blue","#FF6600"],-1,h);