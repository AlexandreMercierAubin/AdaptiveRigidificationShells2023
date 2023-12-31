function mesh3d = shellOBJLoader(filename, attributes, materials, scale, resetMesh, settings)
    %meshLoader(filename, attributes, materials, scale, resetMesh, settings)
    %reads a .mesh file. If possible reads it from a cache mat file
    if nargin < 5
        resetMesh = false;
    end
    
    if nargin < 4
        scale = [1,1,1];
    end

    scalestring = sprintf('%d_',scale);
    tmpPath2 = ['3d/data/cached/', scalestring, materials.cacheName, convertStringsToChars(filename),'.mat'];
    if isfile(tmpPath2) && ~resetMesh
        loadedStruct  = load(tmpPath2);
        mesh3d = loadedStruct.mesh3d;
        if isequal(mesh3d.materials, materials) && isequal(mesh3d.materialIndex, attributes)
            return;
        end
    end
    
    [V,F,UV,TF,N,NF] = readOBJ('3d/data/'+string(filename)+'.obj');
    T = F;
    J = [1:size(T,1)]';
    
%     tetramesh(T,V);
%     vol = volume(V,T);
    if nargin < 4
        scale = [1,1,1];
    end
    V(:,1) = scale(1)*V(:,1);
    V(:,2) = scale(2)*V(:,2);
    V(:,3) = scale(3)*V(:,3);
    mesh3d = Mesh3D(V,T, attributes, materials, T, J,[],[],[],[]);
    save(tmpPath2, 'mesh3d');   
    settings.recomputeCacheAinv = true;
end