clear;
h = 0.01; % time step

% CONFIG
settings = Simulation3DSettings();

% MESHES
cla;
close all;
clear meshes;
rho = 5;
nu = 0.4;      % Poisson ratio: close to 0.5 for rubber
E = 5e4;     % Young's modulus: 0.01e9 approximate for rubber
[ mu, lambda ] = toLame( nu, E );
alpha0 = 0.00001;   % Rayleigh factor on M
alpha1 = 0.000001;  % Rayleigh factor on K
strainUpperBound = 1.3;
strainLowerBound = 0.7;
tMaterial = [TriangleMaterial(rho, mu, lambda, alpha0, alpha1,[1,0,0.5] , strainUpperBound, strainLowerBound)];
thickness = 0.1;
bendingStiffness = materialToBendingStiffness(E,nu,thickness);

radius = 1;
hold on;
th = linspace(0,2*pi,200);
plot3( radius*-(cos(th)-1),radius*sin(th), zeros(size(th)), 'b-');

v1 = [ 0, 0, -0.7 ];
v2 = [ 0, 0,  0.7 ];
t = 1; % initial angle for face size
s = 0.9; % scale on each iteration
num = 50;

for j = 1:num
    v3 = radius*[ -(cos(t)-1), sin(t), 0 ];
    v4 = radius*[ -(cos(-t)-1), sin(-t), 0];
    V = [ v1; v2; v3; v4 ];
    T = [1,2,3;2,1,4];
    J = [1:size(T,1)]';

    baseMesh = Mesh3D(V,T,[],tMaterial,T,J,[],[],[],[]);
    patch('vertices', baseMesh.getPositionFormatted, 'faces', baseMesh.t(:,1:3));
    [curvature, angles, BC, DualEdgeLength, bendingEdgeLength] = computeCurvature(baseMesh, baseMesh.p,1);
    [curvature2, ~] = computeCurvature(baseMesh, baseMesh.p,2);
    [curvature3, ~] = computeCurvature(baseMesh, baseMesh.p,3);

    disp(['t:',num2str(t), ' angle:',num2str(angles), ' curvature - :',num2str(curvature), ' curvature + :',num2str(curvature2)]);

    % shrink edge verts just to make it pretty
    v1 = v1*s;
    v2 = v2*s;
    t = t*s;
end

axis equal;