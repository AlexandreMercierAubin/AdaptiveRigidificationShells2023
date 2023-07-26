clear all;
close all;
h = 0.0025; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesha;
rho = 0.05;
nu = 0.36;      % Poisson ratio: close to 0.5 for rubber
E = 1e0;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
thickness = 0.05;
strainUpperBound = 1.1;
strainLowerBound = 0.90;
kl = 0.01;
ka = 0.5;
bendingAlpha1 = 0.01;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,kl,ka,bendingAlpha1)];

clothDimensions = [2.5,2.5];
precision = 0.04;

mesh = generateCloth([], tMaterial, precision,[1,1,1],clothDimensions);

mesh.setRigidTransform([90,0,0],[0,-0.4,0],true);
mesha = AdaptiveMesh3D(mesh);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 5e1;
rigidificator.ElastificationThreshold = 5e0; 
rigidificator.RigidificationBendThreshold = 5e0;
rigidificator.ElastificationBendThreshold = 5e1;

integrator = LDLBackwardEuler3D();
integrator.useQuicksolveContactFilter = 4;
% integrator = BackwardEuler3D();
% integrator = FullNewton3D();
% integrator.maxIterations = 5;
% integrator.useFullAinv = true;
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0,0 );
% integrator.projectToSPD = true;
% integrator.useFullAinv = true;
% energyModel = StVenantKirchoff3DEnergy();
energyModel = NeoHookean3DEnergy();
% energyModel = CorotationalEnergy();

sphereContactFinder = MovingSphereContactFinder(0.4, [0,-0.3,-0.45], 0.9,h);
sphereContactFinder.plotRatio = 0.95;
sphereContactFinder.cfun = @(t) [0,0,0];
sphereContactFinder.dcdt = @(t) [0,0,0]; 
sphereContactFinder.thetafun = @(t)  .1*(sin(t*h + pi/2) -1)  .* (mod( t*h, pi*2*2 ) < pi*2);
sphereContactFinder.dthetadt = @(t)  .1*cos(t*h + pi/2)     .* (mod( t*h, pi*2*2 ) < pi*2);

planeContactFinder = PlaneContactFinder3D([0,0,1],[0,0,-0.7],0.1);
plotOffset = [0,0,0];
trapezoidCF = ObjectContactDetector("3d/data/trapezoid.obj",[0,0,40], [0.5,0.5,0.25], [1.3,0,-0.7], 0.0, 0.95, plotOffset);
trapezoidCF.faceColor = 'blue';
contactFinder = {sphereContactFinder,planeContactFinder,trapezoidCF};%

settings.MakeVideo = 1;
settings.FramesToRecord = 3500;

% settings.InitialWindowPosition = [0,0,1920,1080];
settings.SceneName = 'clothSphereRotate';
% settings.WriteOBJs = true;
settings.OBJDir = './objs/clothSphereRotate/';
settings.campos=[5,5,2];
settings.PGSiterations = 20;
settings.addBendingEnergy = true;

settings.recomputeCacheAinv = true;
settings.PlotEdotVsCurvatureHists = true;
settings.PlotSkip = plotSkip60FPS(h);
settings.StrainLimitingEnabled = 2;% mexed strain limiting. Use 1 instead if you didn't mex mexStrainLimiting
settings.useGrinspunPlanarEnergy = true;

td = simulate3D(mesha,h,contactFinder, integrator, rigidificator, settings, energyModel);
save("clothSphereRotate_adaptive"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td,"clothSphereRotate_",["adaptive"]);
readTDcsvLog(["clothSphereRotate_default.csv","clothSphereRotate_adaptive.csv"],["blue","#FF6600"],3500,h,true);