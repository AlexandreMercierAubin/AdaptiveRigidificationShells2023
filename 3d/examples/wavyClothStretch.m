%Refer to https://www.tkim.graphics/FEMBW/tkim_sca2020.pdf as a comparison

clear all;
close all;
h = 0.005; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesh3Da;
rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 1e5;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.01;   % Rayleigh factor on M
alpha1 = 0.0001;  % Rayleigh factor on K
thickness = 0.01;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
color = [0.3,0.46,0.8];
ka = 0.5;
kl = 10; 
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka)];


settings.recomputeCacheAinv = true;
dim = [1,1];
sineAlpha = 0.1;
sineFreq = 10;
sinePhase = 0;
length = 2;
mesh3D = generateWavyShell([], tMaterial, 20,length,sineAlpha,sineFreq,sinePhase,[1,1,1]);
mesh3D.setRigidTransform([90,90,0],[0,0,0.5],true);

xPos = mesh3D.p(1:3:end);
yPos = mesh3D.p(2:3:end);
% unpinEnds = (yPos <= max(yPos)-0.2)  & (yPos >= min(yPos)+0.2);
unpinEnds = true(size(yPos)); %disabled unpinned ends
max_I = find(xPos >= max(xPos)-0.1 & unpinEnds);
min_I = find(xPos <= min(xPos)+0.1 & unpinEnds);
mesh3D.pin(sort(max_I));
mesh3D.pin(sort(min_I));

animScripted = SequentialPositionAnimationScripter();
animScripted.dim = 3;

dofs = [min_I*3-2;max_I*3-2];
halfSize= size(dofs,1);
dofsCompress = [min_I*3-1;max_I*3-1];
dofs = [dofs;dofsCompress];

%creates the animated ends
compressFactor = 1;
posInit = mesh3D.p(dofs);
curPos = posInit;
entries = 0.20/h;% hang in the air
start = 0.0/h;
animScripted.frameNumbers = 1:start+entries-1;
posDelta = 0.25;
stackPosChange = [min_I*0-posDelta;max_I*0+posDelta];
compressionDelta = (posInit((halfSize+1):end)*compressFactor - posInit((halfSize+1):end));
counter = 1;
for i = 1:start
    animScripted.dofs{counter} = dofs;
    animScripted.positions{counter} = posInit;
    counter = counter +1;
end
for i = 1:entries
    proportion = i/entries;
    animScripted.dofs{counter} = dofs;
    curPos(1:halfSize) = posInit(1:halfSize)+stackPosChange*proportion;
    curPos((halfSize+1):end) = posInit((halfSize+1):end)+compressionDelta*proportion;
    animScripted.positions{counter} = curPos;
    counter = counter +1;
end

mesh3Da = AdaptiveMesh3D(mesh3D);

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = 0.05;
rigidificator.ElastificationThreshold = 1; 
rigidificator.setBendingThresholdsFromPlanar();
integrator = LDLBackwardEuler3D();
integrator.recordConditionNumber = true;
integrator.Gravity = 0;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.1,0.1 );

energyModel = StVenantKirchoff3DEnergy();
contactFinder = {};

settings.MakeVideo = 1;
settings.SceneName = 'WavyCloth';
settings.OBJDir = './objs/WavyCloth/';
settings.campos=[4,4,5];
settings.camtarget = [0.3,-0.5,0];
settings.camLightPosition = 'left';
settings.addBendingEnergy = true;
settings.useGrinspunPlanarEnergy = true;

simulate3D(mesh3Da,h,contactFinder, integrator, rigidificator, settings, energyModel,animScripted);

compPlot = figure;
ax1 = axes('Parent', compPlot);
hold(ax1,"off");
x = [1:numel(integrator.conditionNumberList)].*h;
semilogy(ax1,x,integrator.conditionNumberList, 'Color','blue');
legend('elastic');
xlabel('Time (seconds)') ;
ylabel('Condition Number Log10') ;