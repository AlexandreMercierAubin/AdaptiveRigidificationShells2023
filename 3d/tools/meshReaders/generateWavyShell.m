function mesh3D = generateWavyShell(attributes, materials, nSides, length, alpha, frequency, phase, meshScale)
    edgeLength = length/nSides;
    layers = ceil(length / edgeLength);
    points = [1:nSides]*edgeLength;

    Vinit = [ points', zeros(nSides,1), zeros(nSides,1) ];
    V = Vinit;
    quad = [1,nSides+1, 2;
            2, nSides+1, nSides+2];
    Tpattern = [];
    for i=0:nSides-2
        curQuad = quad + i;
        Tpattern = [Tpattern;curQuad];
    end

    T=[];
    for i = 1:layers
        zpos = i*edgeLength;
        scale = alpha*sin(frequency*(phase+zpos));
        curV = Vinit + [0,scale,zpos];
        V = vertcat(V,curV);
        T = [T;Tpattern+(i-1)*nSides];
    end
    
    V = V .* meshScale;
    J = [1:size(T,1)]';
    mesh3D = Mesh3D(V,T, attributes, materials, T, J,[],[],[],[]);
%     patch('vertices', V, 'faces', T, 'edgecol', 'k',  'facecol', [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeAlpha', .9 );
%     axis equal;
end

