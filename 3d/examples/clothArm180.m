cla;
clear;
close all;
h = 0.0025; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear mesh3Da;

rho = 1;
nu = 0.37;      % Poisson ratio: close to 0.5 for rubber
E = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = h*0.11;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.1;
strainLowerBound = 0.9;
ka = 10;
kl = 0.02;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,ka,kl,0.001)];

resetMesh = false;
scaling = 1;
mesh3D = shellOBJLoader('armClothLong180',[],tMaterial,[scaling,scaling,scaling],resetMesh,settings);
mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 5e-2;
rigidificator.ElastificationThreshold = 1e-1; 
rigidificator.RigidificationBendThreshold = 5e0;
rigidificator.ElastificationBendThreshold = 7e0;
% integrator = FullNewton3D();
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0.0,0.0);
integrator.useQuicksolveContactFilter = 4;
% energyModel = NeoHookean3DEnergy();
energyModel = StVenantKirchoff3DEnergy();

animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/arm120.mot','./3d/data/AnimationData/armPrettyInv/',[0,0,-90],[1,1,1],[0,0,0],0.5,h);
animContactFinder.plotScale = 0.93;
animContactFinder.interpolateOnFrames = (plotSkip60FPS(h)-1)*3;
% animContactFinder.skipdetection = true;
contactFinder = {animContactFinder};% 

settings.MakeVideo = 1;
settings.SceneName = 'clothArm180';
settings.FramesToRecord = 12.5/h;
settings.PCGiterations = 1;
settings.campos=[10,10,0]*2.5;
settings.camtarget = [0,0,0];
settings.recomputeCacheAinv = true;
settings.StrainLimitingEnabled = 2;
settings.renderer = 'opengl';
settings.OBJDir = './objs/clothArm180Elastic/';
settings.PGSiterations = 10;
settings.PlotSkip = plotSkip60FPS(h);
settings.useGrinspunPlanarEnergy = true;
settings.addBendingEnergy = true;


td = simulate3D({mesh3D},h,contactFinder, integrator, rigidificator, settings, energyModel);
writeTDcsv(td,"clothArm180_",["default"]);
readTDcsvLog(["clothArm180_default.csv","clothArm180_adaptive.csv"],["blue","#FF6600"],-1,h,true);