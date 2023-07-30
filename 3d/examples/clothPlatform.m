clear all;
close all;
scaleFactor = 1;
h = 0.01/scaleFactor; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesh3Da;
rho = 5;
nu = 0.38;      % Poisson ratio: close to 0.5 for rubber
E = 1e3;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.001;   % Rayleigh factor on M
alpha1 = 0.001;  % Rayleigh factor on K
thickness = 0.03;
strainUpperBound = 1.3;
strainLowerBound = 0.7;
bendingAlpha = 0.05;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[0.2,0,0.5], strainUpperBound, strainLowerBound,thickness,1,1,bendingAlpha)];

settings.recomputeCacheAinv = true;
clothDimensions = [2.7,1.8];
resolution = 0.025;
% resolution = 0.06;
mesh3D = generateCloth([], tMaterial, resolution,[1,1,1],clothDimensions);

bottle = meshLoader("bottle", [], tMaterial, [0.1,0.1,0.1],false, settings);
bottle.setRigidTransform([90,0,0],[0,-.1,0.15]);

mesh3D.setRigidTransform([90,0,0],[0,0,0], true);

%creates the animated ends
xPos = mesh3D.p(1:3:end);
min_I = find(xPos <= min(xPos)+0.1);
mesh3D.pin(sort(min_I));

animScripted = SequentialPositionAnimationScripter();
animScripted.dim = 3;
dofs = [min_I*3-2;min_I*3];
posInit = mesh3D.p(dofs);
entries = 400*scaleFactor;% hang in the air
goUp = (3.5/5)*entries;
start = 100*scaleFactor;
posDelta = 0.004/scaleFactor;
for i = 1:entries
    animScripted.dofs{i} = dofs;
    modDelta = posDelta;
    if i < (entries-start)/2 %go up, then down
        zPos = min_I*0+(posDelta/12);
    end
    if i > goUp
        modDelta = posDelta*2;
        zPos = min_I*0 + modDelta;
    end
    stackPosChange = [min_I*0-modDelta;zPos];
    posInit = posInit+stackPosChange;
    animScripted.positions{i} = posInit;
end

%couple steps where nothing happens to zero out pinned dof velocities
counter = entries;
for i = 1:5
    animScripted.dofs{counter} = animScripted.dofs{counter-1};
    animScripted.positions{counter} = animScripted.positions{counter-1};
    counter = counter +1;
end
animScripted.frameNumbers = 1:counter-1;

clothElements = size(mesh3D.t,1);
mesh3D.mergeMesh(bottle);

mesh3Dnr = AdaptiveMesh3D(mesh3D);
mesh3Dnr.renderOffset = [4,0,0];
mesh3Dnr.EnableAutomaticRigidification = false; %This is a way to simulate rigid bodies without the oracle's extra computation
mesh3Dnr.AlwaysRigidElements = [clothElements+1:size(mesh3Dnr.t,1)]';

mesh3Da = AdaptiveMesh3D(mesh3D);
mesh3Da.AlwaysRigidElements = [clothElements+1:size(mesh3Da.t,1)]';

Rt = 5e-3;
Et = 5e-2;
rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = Rt;
rigidificator.ElastificationThreshold = Et; 
rigidificator.setBendingThresholdsFromPlanar();

rigidificator2 = ECurvCloth3DRigidificator();
rigidificator2.RigidificationThreshold = Rt;
rigidificator2.ElastificationThreshold = Et; 
rigidificator2.setBendingThresholdsFromPlanar();

integrator = LDLBackwardEuler3D(); %The copies here are not needed, but let's be safe and make sure they don't have bore effects
integrator.Gravity = -9.8;
integrator.setComplianceAndBaumgarteFromERPandCFM(h, 0.0, 0.0 );

integrator2 = LDLBackwardEuler3D();
integrator2.Gravity = -9.8;
integrator2.setComplianceAndBaumgarteFromERPandCFM(h, 0.0, 0.0 );

energyModel = NeoHookean3DEnergy();

frictionCoeff = 0.0;
plotRatio = 0.98;
platformContactFinder = ObjectContactDetector("3d/data/tableTop.obj",[90,0,0], [1,1,1], [0,0,0], frictionCoeff, plotRatio, [0,0,-0.01]);
platformContactFinder.faceColor = 'blue';
meshMeshCollider = MeshSCD3D(0.3);
contactFinder = {platformContactFinder,meshMeshCollider};

settings.MakeVideo = 1;
settings.FramesToRecord = 30/h;
settings.SceneName = 'clothPlatform';
% settings.WriteOBJs = true;
settings.OBJDir = './objs/clothPlatform/';
settings.StrainLimitingEnabled = false;
settings.campos=[5,5,2];
settings.PGSiterations = 15;
settings.addBendingEnergy = true;
settings.PlotSkip = scaleFactor;
settings.useGrinspunPlanarEnergy = true;

% td = simulate3D({mesh3Dnr,mesh3Da},h,contactFinder, {integrator,integrator2}, {rigidificator,rigidificator2}, settings, energyModel, animScripted );
td = simulate3D({mesh3Da},h,contactFinder, {integrator}, {rigidificator}, settings, energyModel, animScripted );
% save("clothMultiSpheres_"+datestr(now,'mm-dd-yyyy_HH-MM')+".mat", 'td');
% writeTDcsv(td,"clothPlatform_", ["default","adaptive"]);
% readTDcsvLog(["clothPlatform_default.csv","clothPlatform_adaptive.csv"],["blue","#FF6600"],-1,h*1.5,false);  
