function mesh = schur_solve_3D(mesh, h, g)

    % compute right hand side
    rhs = zeros(mesh.N * 3, 1);
    % Gravity
    rhs(3:3:end) = h * mesh.mass(3:3:end) * g;

    % Elastic forces
    rhs = rhs + h * computeForces3D(mesh);

    schur_right = mesh.B * mesh.Minv * rhs;

    F = mesh.B * mesh.p;

    params = zeros(size(mesh.t, 1), 2);
    for n = 1:size(mesh.t, 1)
        params(n, 1:2) = [mesh.lambda, mesh.mu];
    end

    C = d2psi_neohookean_dF2(mesh.t, reshape(F, 9, [])', params);

    for c = 1:9:size(C, 1)
        Cinv(c:c + 8, c:c + 8) = pinv(full(C(c:c + 8, c:c + 8)));
    end

    A = mesh.B * mesh.Minv * mesh.B' + Cinv / h / h;
    
    lambda = A \ schur_right;

    deltav = mesh.Minv * (rhs - mesh.B' * lambda);

    mesh.v = mesh.v + deltav; 
    mesh.p = mesh.p + h * mesh.v;   
    mesh.p(mesh.pinnedDOFs) = mesh.p0(mesh.pinnedDOFs); 
end