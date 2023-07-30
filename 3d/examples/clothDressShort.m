cla;
clear;
close all;
h = 0.0025; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
%bottom
rho = 5;
nu = 0.36;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.0002;   % Rayleigh factor on M
alpha1 = 0.0002;  % Rayleigh factor on K
thickness = 0.015;
ka = 1;
kl = 0.1; 
strainUpperBound = 1.3;
strainLowerBound = 0.7;
bendingAlpha1 = 0.002;
tMaterial1 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,bendingAlpha1)];

%top
rho = 5;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 3e5;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
thickness = 0.015;
ka = 5;
kl = 1; 
tMaterial2 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.2,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,bendingAlpha1)];

%middle
rho = 5;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 1e2;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.015;
ka = 10;
kl = 0.5; 
strainUpperBound = 1.3;
strainLowerBound = 0.7;
tMaterial3 = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.6,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,bendingAlpha1)];
tMaterial= [tMaterial1,tMaterial2,tMaterial3];
resetMesh = true;
mesh3D = shellOBJLoader('dressyellowCoarse',[],tMaterial,[1,1,1],resetMesh,settings);

v = reshape(mesh3D.p,3,[]);
n1 = v( :, mesh3D.t(:,1) );
n2 = v( :, mesh3D.t(:,2) );
n3 = v( :, mesh3D.t(:,3) );
elCenter = (1/3) * (n1+n2+n3);

ind =  elCenter(3:3:end) > 0.05; % this is inverted as they are doubles, not logical

attributes = mesh3D.materialIndex;
attributes(ind) = 3;

ind =  elCenter(3:3:end) > 0.125; % this is inverted as they are doubles, not logical
attributes(ind) = 2;

mesh3D.updateMaterials( attributes, tMaterial );

mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-2;
rigidificator.ElastificationThreshold = 1e-1; 
rigidificator.RigidificationBendThreshold = 20;
rigidificator.ElastificationBendThreshold = 200;
integrator = LDLBackwardEuler3D();
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0,0);
energyModel = NeoHookean3DEnergy();

animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/dance.mot','./3d/data/AnimationData/female-objs/',[-90, 0, 0],[1,1,1],[0, 0, -1],0.9,h);
animContactFinder.plotScale = 0.96;
animContactFinder.interpolateOnFrames = plotSkip60FPS(h);
animContactFinder.animationPositionStart = 17;
animContactFinder.skipPart([1,2,3,4,5,13,14,15,16]) = true;%6,7,8 are the body top, 9,10 are thighs, 11,12 are lower legs 
framePause = ([676,1002, 1206, 1419,1829, 2037,2260, 2539, 2863,3016,3149,3307,3417]-(animContactFinder.interpolateOnFrames*animContactFinder.animationPositionStart)).*(0.0025/h);%evaluated at 0.0025
durations = repmat([floor(0.9/h)],numel(framePause),1);
durations([1]) = floor(0.3/h);
durations([2]) = floor(0.5/h);
durations([5]) = floor(0.7/h);
durations([6]) = floor(0.6/h);
durations([7]) = floor(1/h);
durations([8]) = floor(1.2/h);
animContactFinder.injectPause(framePause,durations);
animContactFinder.stopFrame = (9/h)+sum(durations);
contactFinder = {animContactFinder};

settings.MakeVideo = 1;
settings.SceneName = 'dressShortDefault';
settings.OBJDir = './objs/dressShortDefault/';
settings.FramesToRecord = (10/h)+sum(durations);
settings.PGSiterations = 25;
settings.campos=[-5,-5,0.4];
settings.camtarget = [0,0,0.15];
settings.recomputeCacheAinv = true;
settings.renderer = 'opengl';
settings.camLightPosition = 'headlight';
settings.PlotSkip = plotSkip60FPS(h);
settings.addBendingEnergy = true;
settings.useGrinspunPlanarEnergy =true;

td = simulate3D({mesh3Da},h,contactFinder, integrator, rigidificator, settings, energyModel);
save("clothDressShort_default"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td,"clothDressShort_",["adaptive"]);
readTDcsvLog(["clothDressShort_default.csv","clothDressShort_adaptive.csv"],["blue","#FF6600"],18.7/h,h);