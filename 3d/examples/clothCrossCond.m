%Refer to https://www.tkim.graphics/FEMBW/tkim_sca2020.pdf as a comparison

clear all;
close all;
h = 0.0025; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
clear mesh3Da;
rho = 2;
nu = 0.39;      % Poisson ratio: close to 0.5 for rubber
E = 5e5;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.01;   % Rayleigh factor on M
alpha1 = 0.05;  % Rayleigh factor on K
thickness = 0.25;
strainUpperBound = 1.2;
strainLowerBound = 0.8;
color = [0.3,0.46,0.8];
ka = 125;
kl = 100; 
bendingAlpha1 = 0.00455;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1, [0.3,0.46,0.8], strainUpperBound, strainLowerBound, thickness, kl, ka,bendingAlpha1)];

attributes = [];
materials = tMaterial; 
precision = 0.15;
scale= [1,1,1];
perturbationNormalScaling = 0.001;

sP = 4;
PlusBounds = [-sP,1;
    -1,1;
     -1,sP;
     1,sP;
     1,1;
     sP,1;
     sP,-1;
     1,-1;
     1,-sP;
     -1,-sP;
     -1,-1;
     -sP,-1];

fd = @(p) ddiff(ddiff(ddiff(ddiff(drectangle(p,-sP,sP,-sP,sP),drectangle(p,-sP,-1,-sP,-1)) , drectangle(p,1,sP,1,sP)) ,drectangle(p,-sP,-1,1,sP)), drectangle(p,1,sP,-sP,-1));
[p, T] = distmesh2d(fd, @huniform, precision,[-sP,-sP;sP,sP], PlusBounds);

sumpT1 = abs(p(T(:,1),1) - p(T(:,2),1)) + abs(p(T(:,1),1) - p(T(:,3),1));
isLine1 = sumpT1 <= 1e-8 ;
sumpT2 = abs(p(T(:,1),2) - p(T(:,2),2)) + abs(p(T(:,1),2) - p(T(:,3),2));
isLine2 = abs(sumpT2./3 - p(T(:,1),2)) <= 1e-8 ;

%removes elements that are in the shape of a line
T(isLine1|isLine2, :) = [];

V = [p(:,1),zeros(size(p,1),1),p(:,2)];
V(:,1) = scale(1)*V(:,1);
V(:,2) = scale(2)*V(:,2);
V(:,3) = scale(3)*V(:,3);
J = [1:size(T,1)]';
mesh3D = Mesh3D(V,T, attributes, materials, T, J,[],[],[],[]);

mesh3D.setRigidTransform([90,0.1,0],[0,0,0.5],true);

xPos = mesh3D.p(1:3:end);
yPos = mesh3D.p(2:3:end);

isPinned = false(mesh3D.N,1); %pin ends
distEnds = 0.3;
max_I = find(xPos >= max(xPos)-distEnds);
min_I = find(xPos <= min(xPos)+distEnds);
isPinned(max_I) = true;
isPinned(min_I) = true;
pinVerts = find(isPinned);
mesh3D.pin(pinVerts);

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
entries = 0.5/h;% hang in the air
start = 0.0/h;

posDelta = 3.5;
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

%couple steps where nothing happens to zero out pinned dof velocities
for i = 1:5
    animScripted.dofs{counter} = animScripted.dofs{counter-1};
    animScripted.positions{counter} = animScripted.positions{counter-1};
    counter = counter +1;
end
animScripted.frameNumbers = 1:counter-1;

mesh3D2 = Mesh3D(mesh3D);
mesh3D2.renderOffset = [0,16,0];
mesh3Da = AdaptiveMesh3D(mesh3D);
mesh3Da.renderOffset = [0,-8,0];
mesh3Da2 = AdaptiveMesh3D(mesh3D);
mesh3Da2.renderOffset = [0,8,0];
mesh3Da3 = AdaptiveMesh3D(mesh3D);
mesh3Da3.renderOffset = [0,-16,0];

Et = 5e-3;
Er = 5e-4;

rigidificator = ECurvCloth3DRigidificator();
rigidificator.RigidificationThreshold = Er;
rigidificator.ElastificationThreshold = Et; 
rigidificator.setBendingThresholdsFromPlanar();

rigidificator2 = ECurvCloth3DRigidificator();
rigidificator2.RigidificationThreshold = Er;
rigidificator2.ElastificationThreshold = Et; 
rigidificator2.setBendingThresholdsFromPlanar();

rigidificator3 = ECurvCloth3DRigidificator();
rigidificator3.RigidificationThreshold = Er;
rigidificator3.ElastificationThreshold = Et; 
rigidificator3.setBendingThresholdsFromPlanar();

rigidificator4 = ECurvCloth3DRigidificator();
rigidificator4.RigidificationThreshold = Er;
rigidificator4.ElastificationThreshold = Et; 
rigidificator4.setBendingThresholdsFromPlanar();

rigidificator5 = ECurvCloth3DRigidificator();
rigidificator5.RigidificationThreshold = Er;
rigidificator5.ElastificationThreshold = Et; 
rigidificator5.setBendingThresholdsFromPlanar();

ERP = 0.0;
CFM = 0.0;
iEn = BackwardEuler3D();
iEn.Gravity = -9.8;
iEn.setComplianceAndBaumgarteFromERPandCFM(h, ERP,CFM );
iEn.recordConditionNumber = 1;
iEn.separateQuicksolveGravity = false;

iEs = LDLBackwardEuler3D();
iEs.Gravity = -9.8;
iEs.setComplianceAndBaumgarteFromERPandCFM(h, ERP,CFM );
iEs.recordConditionNumber = 3;
iEs.separateQuicksolveGravity = false;

iAn = BackwardEuler3D();
iAn.Gravity = -9.8;
iAn.setComplianceAndBaumgarteFromERPandCFM(h, ERP,CFM );
iAn.recordConditionNumber = 1;
iAn.separateQuicksolveGravity = false;
% iAn.useFullAinv = true;

iAs = LDLBackwardEuler3D();
iAs.Gravity = -9.8;
iAs.setComplianceAndBaumgarteFromERPandCFM(h, ERP,CFM );
iAs.recordConditionNumber = 3;
iAs.separateQuicksolveGravity = false;
% iAs.useFullAinv = true;

iAsJ = BackwardEuler3D();
iAsJ.Gravity = -9.8;
iAsJ.setComplianceAndBaumgarteFromERPandCFM(h, ERP,CFM );
iAsJ.recordConditionNumber = 2; %2 is jacobi, 5 is block jacobi
iAsJ.separateQuicksolveGravity = false;

energyModel = StVenantKirchoff3DEnergy();
contactFinder = {};

settings.MakeVideo = 1;
settings.FramesToRecord = 3/h;
settings.PlotSkip = plotSkip60FPS(h);
settings.SceneName = 'clothStretchCond';
settings.WriteOBJs = true;
settings.OBJDir = './objs/clothStretchCond/';
settings.campos=[7,10,3]*4;
settings.camLightPosition = 'left';
settings.addBendingEnergy = true;
settings.useGrinspunPlanarEnergy = true;
settings.StrainLimitingEnabled = true;
settings.recomputeCacheAinv = true;

simulate3D({mesh3D,mesh3D2,mesh3Da,mesh3Da2,mesh3Da3},h,contactFinder,  {iEn,iEs,iAn,iAs,iAsJ}, {rigidificator,rigidificator2,rigidificator3,rigidificator4,rigidificator5}, settings, energyModel,animScripted);

fontName = 'Linux Biolinum O';
compPlot = figure('Renderer', 'painters', 'Units','Inches', 'Position', [0,0,3.3125,2]);
fontSize = 8;
set(gcf,'color','w');
ax1 = axes('Parent', compPlot);
hold(ax1,"off");
range = 1:numel(iEn.conditionNumberList);
x = [1:numel(iEn.conditionNumberList(range))].*h;
semilogy(ax1,x,iEn.conditionNumberList(range), 'Color','#9999FF');
hold(ax1,"on");
piEs = semilogy(ax1,x,iEs.conditionNumberList(range), 'Color','#152299', "LineWidth",1);
% piEs.LineStyle = "-";
semilogy(ax1,x,iAn.conditionNumberList(range), 'Color',"#FF6600");
piAn = semilogy(ax1,x,iAs.conditionNumberList(range), 'Color',"#EE0000", "LineWidth",1);
semilogy(ax1,x,iAsJ.conditionNumberList(range), 'Color','#710193');
% piAn.LineStyle = "-";
legend('Elastic','Elastic+ALS','Adaptive', 'Adaptive+ALS', 'Adaptive+Jacobi');
xlabel('Time (seconds)') ;
ylabel('Condition Number Log10') ;
ylim([1,max(iAn.conditionNumberList(range))]);
xlim([0,max([1:numel(iEn.conditionNumberList(range))].*h)]);
set(gca, 'YTick', [10.^0 10.^2 10.^4 10.^6]);
handle = gca;
handle.Color = [0.95, 0.95, 0.95];
handle.TickLength = [0,0];
handle.YGrid = 'on';
handle.GridColor = [1 1 1];
handle.XColor = [0 0 0];
set(gca, 'fontSize', fontSize, 'fontName', fontName);