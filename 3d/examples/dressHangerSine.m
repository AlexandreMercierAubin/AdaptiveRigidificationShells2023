cla;
clear;
close all;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

clear mesh3Da;

rho = 5;
nu = 0.35;      % Poisson ratio: close to 0.5 for rubber
E = 2e3;     % Young's modulus: 0.01e9 approximate for rubber
[mu, lambda] = toLame(nu, E);
alpha0 = 0.0005;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.01;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
ka = 20;
kl = 0.5;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,kl,ka)];

resetMesh = true;
mesh3D = shellOBJLoader('dressyellowCoarse',[],tMaterial,[1,1,1],resetMesh,settings);
mesh3D.setRigidTransform([0,0,90],[0,0,1.2],true);

% baseMesh.pin([1,3,4]);

mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-5;
rigidificator.ElastificationThreshold = 1e-4; 
rigidificator.setBendingThresholdsFromPlanar();
integrator = LDLBackwardEuler3D();
integrator.useQuicksolveContactFilter = 4;
integrator.setComplianceAndBaumgarteFromERPandCFM(h,0,0.0);
energyModel = StVenantKirchoff3DEnergy();

mannequinPosition = [0,0,1.5];
plotRatio= 1;
frictionCoeff=0.7;
animContactFinder = AnimatedObjectContactDetector('./3d/data/AnimationData/platformSine.mot','./3d/data/AnimationData/hanger/',[0,0,0],[1.5,0.4,0.2],mannequinPosition,frictionCoeff,h);
animContactFinder.plotScale = plotRatio;
animContactFinder.interpolateOnFrames = 0.04/h;

meshMeshContactFinder = MeshSCD3D(0.0,false, []);
contactFinder = {animContactFinder,meshMeshContactFinder};

settings.MakeVideo = 1;
settings.SceneName = 'hangerScene';
settings.PCGiterations = 1;
settings.PGSiterations = 30;
settings.campos=[5,5,1];
settings.camtarget = [0,0,1];
settings.recomputeCacheAinv = true;
settings.useGrinspunPlanarEnergy = true;
settings.useGrinspunPlanarEnergy = true;
settings.addBendingEnergy = true;
settings.PlotSkip = plotSkip60FPS(h);

td = simulate3D({mesh3Da},h,contactFinder, integrator, rigidificator, settings, energyModel);

compPlot = figure;
ax1 = axes('Parent', compPlot);
hold(ax1,"off");
x = [1:numel(integrator.conditionNumberList)].*h;
semilogy(ax1,x,integrator.conditionNumberList, 'Color','cyan');
legend('elastic');
xlabel('Time (seconds)') ;
ylabel('Condition Number') ;