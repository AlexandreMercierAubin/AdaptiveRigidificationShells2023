%The difference between 1 and 2 is that this one is mostly edge aligned
%with the bending direction
clear;
cla;
clear mesha;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
thickness = 0.1;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
ka = 1999;
kl = 1999; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,0.0001)];

resolutionHD = 0.017;
resolutionMED = 0.03;
resolutionLOW = 0.045;
direction = 2;
meshLOW = generateCloth([], tMaterial, resolutionLOW,[1,1,1],[3,1],0,direction);
meshMED = generateCloth([], tMaterial, resolutionMED,[1,1,1],[3,1],0,direction);
meshHD = generateCloth([], tMaterial, resolutionHD,[1,1,1],[3,1],0,direction);

settings.recomputeCacheAinv = true;

% pinning tris
xPos = meshLOW.p(1:3:end);
max_I = find(xPos >= max(xPos)-0.5);
meshLOW.pin(sort(max_I));
meshLOW.setRigidTransform([0,90,-90],[0,0,0]);
meshLOW.renderOffset=[-1.1,0,0];
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
meshHD.renderOffset=[1.1,0,0];
meshaHD = AdaptiveMesh3D(meshHD);

Et = 999999;
Rt = 999999;
Eb = 5e-2;
Rb = 5e-3;
bendtype = 3;%1 curvature, 2 dihedral, 3 curvature with rest dual

rigidificationOn = true;
meshaLOW.EnableAutomaticRigidification = rigidificationOn;
meshaMED.EnableAutomaticRigidification = rigidificationOn;
meshaHD.EnableAutomaticRigidification = rigidificationOn;

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = Rt;
rigidificator.ElastificationThreshold = Et; 
rigidificator.ElastificationBendThreshold = Eb;
rigidificator.RigidificationBendThreshold = Rb;
rigidificator.PreventPinnedRigidification = true;
rigidificator.bendType = bendtype;

rigidificator2 = ECurvCloth3DRigidificator();
rigidificator2.RigidificationThreshold = Rt;
rigidificator2.ElastificationThreshold = Et;
rigidificator2.ElastificationBendThreshold = Eb;
rigidificator2.RigidificationBendThreshold = Rb;
rigidificator2.PreventPinnedRigidification = true;
rigidificator2.bendType = bendtype;

rigidificator3 = ECurvCloth3DRigidificator();
rigidificator3.RigidificationThreshold = Rt;
rigidificator3.ElastificationThreshold = Et; 
rigidificator3.ElastificationBendThreshold = Eb;
rigidificator3.RigidificationBendThreshold = Rb;
rigidificator3.PreventPinnedRigidification = true;
rigidificator3.bendType = bendtype;

integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;

integrator2 = LDLBackwardEuler3D();
integrator2.Gravity = -9.8;

integrator3 = LDLBackwardEuler3D();
integrator3.Gravity = -9.8;

energyModel = StVenantKirchoff3DEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.MakeVideo = 1;
settings.FramesToRecord = 30/h;
settings.SceneName = 'clothResTest2';
settings.WriteOBJs = true;
settings.OBJDir = './objs/clothResTest2/';
settings.addBendingEnergy = 1; %Note: the bending energy is very slow to
% compute with matlab. Using a mex version would speed it up.
settings.campos=[9,9,2];
settings.addShellNormalDeformation = true;

td = simulate3D({meshaLOW,meshaMED,meshaHD},h,contactFinder, {integrator,integrator2,integrator3}, {rigidificator,rigidificator2,rigidificator3}, settings, energyModel);
writeTDcsv(td,"clothRes",["low","med","hd"]);
writematrix(rigidificator.curvatures,"curvaturesLOW");
writematrix(rigidificator2.curvatures,"curvaturesMED");
writematrix(rigidificator3.curvatures,"curvaturesHD");