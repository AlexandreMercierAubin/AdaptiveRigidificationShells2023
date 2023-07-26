function mesh = custom_idea_3D(mesh, h, g)

    countLambda = size(mesh.B, 1);

    % compute right hand side
    rhs = zeros(mesh.N * 3 + countLambda, 1);
    % Gravity
    rhs(3:3:(mesh.N * 3)) = h * mesh.mass(3:3:end) * g;

    % Elastic forces
    rhs(1:(mesh.N * 3)) = rhs(1:(mesh.N * 3)) + h * computeForces3D(mesh);


    F = mesh.B * mesh.p;

    params = zeros(size(mesh.t, 1), 2);
    for n = 1:size(mesh.t, 1)
        params(n, 1:2) = [mesh.lambda, mesh.mu];
    end


    C = d2psi_neohookean_dF2(mesh.t, reshape(F, 9, [])', params);
    Acell = { mesh.M, mesh.B'; -h * h * C * mesh.B, eye(countLambda)};

    A = cell2mat(Acell);
    
    rhs(1:(mesh.N * 3)) = rhs(1:(mesh.N * 3)) - h * h * mesh.B' * C * mesh.B * mesh.v;

    ii = cat(2, mesh.unpinnedDOFs, 3 * mesh.N + 1:3 * mesh.N + countLambda);
    
    deltav(ii) = A(ii, ii) \ rhs(ii);

    mesh.v = mesh.v + deltav(1:(3 * mesh.N))'; 
    mesh.p = mesh.p + h * mesh.v;   
    mesh.p(mesh.pinnedDOFs) = mesh.p0(mesh.pinnedDOFs); 
end