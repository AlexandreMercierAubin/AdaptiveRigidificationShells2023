clear;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesh3Da;
rho = 10;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 5e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.01;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness,5,1,0.001)];

resolution = 0.1;
mesh3D = generateCloth([], tMaterial, resolution,[1,1,1],[3,1],0.00001);

% pinning tris
zPos = mesh3D.p(3:3:end);
xPos = mesh3D.p(1:3:end);
max_I = find(zPos >= max(zPos)-0.01 & (xPos == max(xPos) | xPos == min(xPos)));
mesh3D.pin(sort(max_I));

%creates the animated ends
animScripted = SequentialPositionAnimationScripter();
animScripted.dim = 3;
dofs = max_I*3-2;
multiplyer = 1/100/h;
posInit = mesh3D.p(dofs);
entries = 50*multiplyer;% hang in the air
start = 10*multiplyer;
animScripted.frameNumbers = start:start+entries-1;
posDelta = 0.006*(1/multiplyer);
for i = 1:start
    animScripted.dofs{i} = dofs;
    animScripted.positions{i} = posInit;
end
for i = start:entries
    animScripted.dofs{i} = dofs;
    stackPosChange = [posDelta;-posDelta];
    posInit = posInit+stackPosChange;
    animScripted.positions{i} = posInit;
end

mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 1e-4;
rigidificator.ElastificationThreshold = 5e-3; 
integrator = LDLBackwardEuler3D();
integrator.recordConditionNumber = true;
integrator.Gravity = -9.8;

energyModel = StVenantKirchoff3DEnergy();

NullContactFinder = NullContactFinder(3);
contactFinder = {NullContactFinder};

settings.StrainLimitingEnabled = false;
settings.MakeVideo = 1;
settings.SceneName = 'hangingClothStillTwoPin';
settings.addShellNormalDeformation = 1;
settings.addBendingEnergy = true;
settings.campos=[5,5,0];
settings.recomputeCacheAinv = true;
settings.PlotEdotVsCurvatureHists = true;
settings.RigidificationEnabled = false;
settings.elementLineColor = 'black';
settings.useGrinspunPlanarEnergy = true;

simulate3D(mesh3Da,h,contactFinder, integrator, rigidificator, settings, energyModel,animScripted);

compPlot = figure;
ax1 = axes('Parent', compPlot);
hold on;
x = [1:numel(integrator.conditionNumberList)].*h;
semilogy(ax1,x,integrator.conditionNumberList, 'Color','cyan');
legend('elastic');
xlabel('Time (seconds)') ;
ylabel('Condition Number') ;