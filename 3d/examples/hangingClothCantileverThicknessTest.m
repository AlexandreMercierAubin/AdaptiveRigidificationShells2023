clear;
cla;
clear mesha;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
strainUpperBound = 1.3;
strainLowerBound = 0.7;
ka = 5;
kl = 10; 
tMaterialH = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, 0.2, kl, ka)];
tMaterialM = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, 0.1, kl, ka)];
tMaterialL = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, 0.05, kl, ka)];

resolutionMED = 0.05;
meshLOW = generateCloth([], tMaterialL, resolutionMED,[1,1,1]);
meshMED = generateCloth([], tMaterialM, resolutionMED,[1,1,1]);
meshHD = generateCloth([], tMaterialH, resolutionMED,[1,1,1]);
settings.recomputeCacheAinv = true;

% pinning tris
xPos = meshLOW.p(1:3:end);
max_I = find(xPos >= max(xPos)-0.5);
meshLOW.pin(sort(max_I));
meshLOW.setRigidTransform([0,90,-90],[0,0,0]);
meshLOW.renderOffset=[1.1,0,0];
meshaLOW = AdaptiveMesh3D(meshLOW);


xPos = meshMED.p(1:3:end);
max_I = find(xPos >= max(xPos)-0.5);
meshMED.pin(sort(max_I));
meshMED.setRigidTransform([0,90,-90],[0,0,0]);
meshaMED = AdaptiveMesh3D(meshMED);

xPos = meshHD.p(1:3:end);
max_I = find(xPos >= max(xPos)-0.5);
meshHD.pin(sort(max_I));
meshHD.setRigidTransform([0,90,-90],[0,0,0]);
meshHD.renderOffset=[-1.1,0,0];
meshaHD = AdaptiveMesh3D(meshHD);

Et = 5e-2;
Rt = 5e-3;
rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = Rt;
rigidificator.ElastificationThreshold = Et; 
rigidificator.setBendingThresholdsFromPlanar();
rigidificator.PreventPinnedRigidification = true;

rigidificator2 = ECurvCloth3DRigidificator();
rigidificator2.RigidificationThreshold = Rt;
rigidificator2.ElastificationThreshold = Et; 
rigidificator2.setBendingThresholdsFromPlanar();
rigidificator2.PreventPinnedRigidification = true;

rigidificator3 = ECurvCloth3DRigidificator();
rigidificator3.RigidificationThreshold = Rt;
rigidificator3.ElastificationThreshold = Et; 
rigidificator3.setBendingThresholdsFromPlanar();
rigidificator3.PreventPinnedRigidification = true;


integrator = LDLBackwardEuler3D();
% integrator.maxIterations = 5;
integrator.Gravity = -9.8;

integrator2 = LDLBackwardEuler3D();
% integrator2.maxIterations = 5;
integrator2.Gravity = -9.8;

integrator3 = LDLBackwardEuler3D();
% integrator3.maxIterations = 5;
integrator3.Gravity = -9.8;

energyModel = StVenantKirchoff3DEnergy();
% energyModel = NeoHookean3DEnergy();
% energyModel = CorotationalEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.PlotEdotVsCurvatureHists = true;
settings.MakeVideo = 1;
settings.FramesToRecord = 30/h;
% settings.PlotEDotHist = 1;
% settings.InitialWindowPosition = [0,0,1920,1080];
settings.SceneName = 'clothThickTest';
settings.WriteOBJs = true;
settings.OBJDir = './objs/clothThickTest/';
% settings.StrainLimitingEnabled = true;
settings.addBendingEnergy = 1; %Note: the bending energy is very slow to
% compute with matlab. Using a mex version would speed it up.
settings.campos=[9,9,2];
% settings.RigidificationEnabled = false;
% settings.PlotDihedralRateHist = true;
settings.useGrinspunPlanarEnergy = true;
% settings.addShellNormalDeformation = true;

td = simulate3D({meshaLOW,meshaMED,meshaHD},h,contactFinder, {integrator,integrator2,integrator3}, {rigidificator,rigidificator2,rigidificator3}, settings, energyModel);
writeTDcsv(td,"clothThick",["low","med","hd"]);
