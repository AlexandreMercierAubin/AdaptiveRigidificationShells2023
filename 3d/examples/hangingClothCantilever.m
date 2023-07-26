clear;
cla;
clear mesha;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES

rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 3e5;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.0001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.05;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
ka = 10;
kl = 10; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka)];


% resolution = 0.05;
resolution = 0.2;
direction = 2;
mesh = generateCloth([], tMaterial, resolution,[1,1,1],[3,1],0,direction);
settings.recomputeCacheAinv = true;

% pinning tris
xPos = mesh.p(1:3:end);
max_I = find(xPos >= max(xPos)-0.5);
mesh.pin(sort(max_I));

mesh.setRigidTransform([0,90,-90],[0,0,0]);
mesha = AdaptiveMesh3D(mesh);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 5e-4;
rigidificator.ElastificationThreshold = 5e-3; 
rigidificator.setBendingThresholdsFromPlanar(3);
% rigidificator.PreventPinnedRigidification = true;%might want to turn this
% one if things lock too quickly as pinned + rigidified dofs are completely
% removed from the system

integrator = LDLBackwardEuler3D();
integrator.Gravity = -9.8;
energyModel = StVenantKirchoff3DEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.PlotEdotVsCurvatureHists = true;
settings.MakeVideo = 1;
settings.SceneName = 'hangingClothCantilever';
settings.addBendingEnergy = 1; %Note: the bending energy is very slow to
% compute with matlab. Using a mex version would speed it up.
settings.campos=[9,9,2];
settings.useGrinspunPlanarEnergy = true;

simulate3D(mesha,h,contactFinder, integrator, rigidificator, settings, energyModel);

figure;
x = [1:numel(integrator.conditionNumberList)].*h;
plot(x,integrator.conditionNumberList);
xlabel('Time (seconds)') ;
ylabel('Condition Number') ;
