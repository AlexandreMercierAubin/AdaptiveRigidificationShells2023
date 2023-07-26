function mesh3D = generateShellTube(attributes, materials, nSides, length, initialRadius, alpha, frequency, phase, meshScale)
%     alpha = 0.1;
%     frequency = 0.2;
    dtheta = 2*pi / nSides;
    theta = 0:dtheta:dtheta*(nSides-1);
    
    s = sin(theta);
    c = cos(theta);
    circleV = [ initialRadius*c', initialRadius*s', zeros(nSides,1) ];
    edgeLength = norm(circleV(1,:) - circleV(2,:));
    layers = ceil(length / edgeLength);

    scale = 1 + alpha*sin(frequency*(phase));
    V = circleV.*[scale,scale,1] + [0,0,0];

    quad = [1,nSides+1, 2;
            2, nSides+1, nSides+2];
    Tpattern = [];
    for i=0:nSides-2
        curQuad = quad + i;
        Tpattern = [Tpattern;curQuad];
    end

    %close the loop
    Tpattern = [Tpattern;
        nSides,2*nSides,1;
        1, nSides*2,nSides+1];
    
    T=[];
    for i = 1:layers
        zpos = i*edgeLength;
        scale = 1 + alpha*sin(frequency*(phase+zpos));
        curV = circleV.*[scale,scale,1] + [0,0,zpos];
        V = vertcat(V,curV);
        T = [T;Tpattern+(i-1)*nSides];
    end
    
    V = V .* meshScale;
    J = [1:size(T,1)]';
    mesh3D = Mesh3D(V,T, attributes, materials, T, J,[],[],[],[]);
%     patch('vertices', V, 'faces', T, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeAlpha', .9 );
%     axis equal;
end

