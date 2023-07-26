clear all;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesha;
rho = 10;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 2e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.05;
strainLowerBound = 0.95;
kl = 1;
ka = 1;
bendingAlpha1 = 0.001;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,kl,ka,bendingAlpha1)];

clothDimensions = [3,3];
% clothDimensions = [3,1.5];
precision = 0.1;%0.035 for final
mesh = generateCloth([], tMaterial, precision,[1,1,1],clothDimensions);

mesh.setRigidTransform([90,0,0],[0,-0.3,0.5]);
mesha = AdaptiveMesh3D(mesh);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 1e-3; 
rigidificator.setBendingThresholdsFromPlanar(3);

integrator = LDLBackwardEuler3D();

integrator.Gravity = -9.8;
integrator.useQuicksolveContactFilter = 4;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0,0 );
energyModel = NeoHookean3DEnergy();

sphereContactFinder = SphereContactFinder(0.7, [0,-0.5,-0.3], 0.8);
sphereContactFinder.plotRatio = 0.95;

contactFinder = {sphereContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 300/h;

settings.SceneName = 'clothSphere';
settings.OBJDir = './objs/clothSphere/';
settings.campos=[5,5,2];
settings.PGSiterations = 20;
settings.addBendingEnergy = true;

settings.recomputeCacheAinv = true;
settings.PlotEdotVsCurvatureHists = true;

settings.PlotSkip = plotSkip60FPS(h);
settings.StrainLimitingEnabled = 0;
settings.useGrinspunPlanarEnergy = true;

td = simulate3D(mesha,h,contactFinder, integrator, rigidificator, settings, energyModel);