clear;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear meshes;
rho = 2;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.01;  % Rayleigh factor on K
thickness = 0.005;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
ka = 0.01;
kl = 0.01;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,kl,ka,0.001)];

baseMesh = generateCloth([], tMaterial, 0.05);

% pinning tris
zPos = baseMesh.p(3:3:end);
max_I = find(zPos >= max(zPos)-0.01);
baseMesh.pin(sort(max_I));

meshes = AdaptiveMesh3D(baseMesh);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 1e-3; 
integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;
integrator.separateQuicksolveGravity = false;

energyModel = StVenantKirchoff3DEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.StrainLimitingEnabled = false;
settings.MakeVideo = 1;
settings.SceneName = 'hangingClothStill';
settings.WriteOBJs = true;
settings.OBJDir = './objs/hangingClothStill/';
settings.addBendingEnergy = true;
settings.campos=[5,5,0];
settings.PlotEdotVsCurvatureHists = true;
settings.useGrinspunPlanarEnergy = true;
settings.recomputeCacheAinv = true;

simulate3D(meshes,h,contactFinder, integrator, rigidificator, settings, energyModel);