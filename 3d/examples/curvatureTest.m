
clear

figure(2)
clf

th = linspace(0,2*pi,300);
subplot(1,2,1);
radius = 1;
plot3( radius*-(cos(th)-1),radius*sin(th), zeros(size(th)), 'b-');
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
hold on;

v1 = [ 0, 0, -0.7 ];
v2 = [ 0, 0,  0.7 ];
t = 1; % initial angle for face size
s = 0.9; % scale on each iteration
num = 50;
kappa1 = zeros(num,1);
kappa2 = zeros(num,1);
kappa3 = zeros(num,1);
kappa4 = zeros(num,1);

for j = 1:num
    v3 = radius*[ -(cos(t)-1), sin(t), 0 ];
    v4 = radius*[ -(cos(-t)-1), sin(-t), 0];
    V = [ v1; v2; v3; v4 ];
    T = [1,2,3;2,1,4];
    tm = trimesh( T, V(:,1), V(:,2), V(:,3) );
    tm.FaceAlpha = 0;
    tm.EdgeColor = [0,0,0];

    N = faceNormals( T, V );
    n1 = N(1,:);
    n2 = N(2,:);
    % make a local frame at the edge
    e1 = makeUnitLength(v2 - v1);
    e2 = n1;
    e3 = makeUnitLength(cross( e1, e2 ));
    E = [e1', e2', e3'];
    % find coords of N2 in this frame
    c = E'*n2';
    theta1 = atan2( c(3), c(2) );   
    % the unsigned version of the angle to debug
    % theta2 = acos( dot(n1,n2) ); 

    BC = barycenters(T,V);
%     BC = barycenter(V,T);

    % should dual edge run along the surface?
    % should it be projected to be perp to normal?
    % all this might be irrelevant for small triangles
    dualEdge = BC(2,:) - BC(1,:);
    dualEdgeLen = norm(dualEdge);

    % here is the geodesic distance
    % could do this sum to the midpoint of the edge
    % to make an easy version that works in general.
    dualEdgeGeoLen = norm(BC(2,:)) + norm(BC(1,:));

    del(j) = dualEdgeLen;
    kappa1(j) = 4 * tan(theta1/2) / dualEdgeGeoLen;
    kappa2(j) = 2* theta1 / dualEdgeGeoLen;
    kappa3(j) = 4 * tan(theta1/2) / dualEdgeLen;
    kappa4(j) = 2* theta1 / dualEdgeLen
    
    % shrink edge verts just to make it pretty
    v1 = v1*s;
    v2 = v2*s;
    t = t*s;
end
subplot(1,2,2);
plot( 1:num, kappa1,'b');
hold on;
plot( 1:num, kappa2,'r');
plot( 1:num, kappa3,'g');
plot( 1:num, kappa4,'k');
legend(["1","2","3","4"]);
yline( 1/radius, 'k:')
ylabel('kappa estimate')'
xlabel('smaller and smaller triangles')
% look at the region around the true curvature
ylim([1/radius/1.2,1.2*1/radius]);

function n = faceNormals( T, V )
    e1 = V(T(:,2),:) - V(T(:,1)); 
    e2 = V(T(:,3),:) - V(T(:,1)); 
    n = cross( e1, e2 );
    n = makeUnitLength(n);
end

function v = makeUnitLength( v )
    % makes vectors stored in different rows unit length
    v = v ./ sqrt( sum ( v.^2, 2 ) );
end

function BC = barycenters( T, V )
    BC = zeros( size(T) );    
    for i = 1:size(T,1)
        BC(i,:) = sum( V(T(i,:),:,1));
    end
end