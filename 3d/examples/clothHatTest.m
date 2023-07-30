clear;
close all;
cla;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
rho = 1;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 1e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0005;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.02;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
color = [0.3,0.46,0.8];
ka = 1;
kl = 5; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka)];

%TOP
rho = 2;
ka2 = 10;
kl2 = 50; 
thickness = 0.1;
nu2 = 0.4;      % Poisson ratio: close to 0.5 for rubber
E2 = 5e4;
[ mu2, lambda2 ] = toLame( nu2, E2 );
alpha0 = 0.005;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
tMaterial2 = [TriangleMaterial(rho, mu2, lambda2, alpha0, alpha1, [0.3,0.2,0.8], strainUpperBound, strainLowerBound, thickness, kl2, ka2)];
tMaterial= [tMaterial,tMaterial2];


resetMesh = true;
mesh3D = shellOBJLoader('hathd',[],tMaterial,[0.5,0.5,0.5],resetMesh,settings);


v = reshape(mesh3D.p,3,[]);
n1 = v( :, mesh3D.t(:,1) );
n2 = v( :, mesh3D.t(:,2) );
n3 = v( :, mesh3D.t(:,3) );
elCenter = (1/3) * (n1+n2+n3);
% The .24 here is because it doesn't quite sit on the y axis... small shift
% in the z direction.
dfromyaxis = sqrt(sum((elCenter([1,3],:)-[0;0.05]).^2, 1));
dfromyaxis2 = sqrt(sum((elCenter([1,3],:)-[0;-0.1]).^2, 1));
rthresh = 0.3;
rthresh2 = 0.34;
% find the bits away from the core, and avoid the boxes on the side of teh
% rocket or other features that are too far off the centerline
ind =  (dfromyaxis > rthresh  & abs(elCenter(2,:))<1) &(dfromyaxis2 > rthresh2 & abs(elCenter(2,:))<1); % this is inverted as they are doubles, not logical

attributes = mesh3D.materialIndex;
attributes(ind) = 2;
mesh3D.updateMaterials( attributes, tMaterial );

% generateCloth([], tMaterial, 0.1);
mesh3D.setRigidTransform([70,30,0],[0,0,0.4],true);
mesh3D.renderOffset = [2,0,0];
% pinning tris
mesh3Da = AdaptiveMesh3D(mesh3D);
mesh3Da.renderOffset = [0,0,0];

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-3;
rigidificator.ElastificationThreshold = 1e-2; 
rigidificator.setBendingThresholdsFromPlanar();
rigidificator.bendType = 1;

integrator4 = LDLBackwardEuler3D();
integrator4.Gravity = -9.8;
integrator4.setComplianceAndBaumgarteFromERPandCFM(h, 0.0, 0.0);

energyModel = StVenantKirchoff3DEnergy();

planeContactFinder = PlaneContactFinder3D([0,0,1], [0,0,-0.2], 0.0);
contactFinder = {planeContactFinder};

settings.FramesToRecord = 4.2/h; %time in seconds scaled by h
settings.SceneName = 'hat';
settings.OBJDir = './objs/hat/';
settings.campos=[8,3,1.75];
settings.addBendingEnergy = 1;

settings.recomputeCacheAinv = true;
settings.PGSiterations = 15;
settings.useGrinspunPlanarEnergy = true;
settings.PCGiterations = 2;

td = simulate3D({mesh3Da},h,contactFinder, {integrator4}, rigidificator, settings, energyModel);%
save("hat_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
writeTDcsv(td, "hat", ["_adaptiveFilter","_adaptiveCG","_adaptiveOG","_default"]);
% readTDcsv(["hat_adaptive.csv","hat_default.csv"],-1,(h/10)/(1/60));
