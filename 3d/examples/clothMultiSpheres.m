clear all;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesha;
rho = 10;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.005;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
bendingAlpha1 = 0.05;
thickness = 0.01;
strainUpperBound = 1.1;
strainLowerBound = 0.9;
kl = 0.3;
ka = 0.5;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,kl,ka,bendingAlpha1)];

% clothDimensions = [3,3];
clothDimensions = [3,1.5];
precision = 0.025;%0.025 for final
mesh = generateCloth([], tMaterial, precision,[1,1,1],clothDimensions);

mesh.setRigidTransform([90,0,0],[0,-0.3,0.5]);
% mesh.renderOffset = [4,0,0];
mesha = AdaptiveMesh3D(mesh);
mesha.renderOffset = [0,0,0];

rigidificator = ECurvCloth3DRigidificator();
%standard
% rigidificator.RigidificationThreshold = 1e-2;
% rigidificator.ElastificationThreshold = 1e-1; 
% rigidificator.RigidificationBendThreshold = 1e-1;
% rigidificator.ElastificationBendThreshold = 1e-0;

%SLimit
rigidificator.RigidificationThreshold = 7e-2;
rigidificator.ElastificationThreshold = 7e-1; 
rigidificator.RigidificationBendThreshold = 5e-1;
rigidificator.ElastificationBendThreshold = 2e-0;


integrator = LDLBackwardEuler3D();
% integrator = BackwardEuler3D();
% integrator = FullNewton3D();
% integrator.maxIterations = 5;
% integrator.useFullAinv = true;
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.0,0.0 );
% integrator.projectToSPD = true;
% integrator.useFullAinv = true;
% energyModel = StVenantKirchoff3DEnergy();
energyModel = NeoHookean3DEnergy();
integrator.useQuicksolveContactFilter = 4;
% energyModel = CorotationalEnergy();

sphereContactFinder = SphereContactFinder(0.7, [0,-1,-0.3], 0.8);
sphereContactFinder.plotRatio = 0.95;
sphereContactFinder2 = SphereContactFinder(0.4, [-0.5,0,-0.3], 0.8);
sphereContactFinder2.plotRatio = 0.95;
sphereContactFinder3 = SphereContactFinder(0.3, [0.5,0.1,-0.3], 0.8);
sphereContactFinder3.plotRatio = 0.95;
sphereContactFinder4 = SphereContactFinder(0.2, [1,-0.1,-0.3], 0.8);
sphereContactFinder4.plotRatio = 0.94;
contactFinder = {sphereContactFinder,sphereContactFinder2,sphereContactFinder3,sphereContactFinder4};

settings.MakeVideo = 1;
settings.FramesToRecord = 11/h;

% settings.InitialWindowPosition = [0,0,1920,1080];
settings.SceneName = 'clothMultiSphere';
% settings.WriteOBJs = true;
settings.OBJDir = './objs/clothMultiSphere/';
settings.campos=[5,5,2];
settings.PGSiterations = 10;
% settings.quicksolveSimulation = true;
settings.addBendingEnergy = true;

% settings.RigidificationEnabled = false;
settings.recomputeCacheAinv = true;
settings.PlotEdotVsCurvatureHists = true;
% settings.PlotEDotHist = 1;
settings.PlotSkip = plotSkip60FPS(h);
settings.StrainLimitingEnabled = 1;% mexed strain limiting. Use 1 instead if you didn't mex mexStrainLimiting
settings.useGrinspunPlanarEnergy = true;

td = simulate3D({mesh,mesha},h,contactFinder, integrator, rigidificator, settings, energyModel);
save("clothMultiSpheres_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td, "clothMultiSpheres_", ["default","adaptive"]);
readTDcsvLog(["clothMultiSpheres_default.csv","clothMultiSpheres_adaptive.csv"],["blue","#FF6600"],700,h,true);