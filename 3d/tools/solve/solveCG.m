function x = solveCG(iterations, A, b, initial_x, tol)
    if nargin < 5
        tol = 1e-30;
    end    
    x = initial_x;
    r = b - A * x;
    d = r;
    
    for i = 1:iterations
        alpha = (r' * r) / (d' * A * d);
        x = x + alpha * d;
        
        oldr = r;
        r = r - alpha * A * d;
        if norm(r) < tol
              break;
        end
        beta = (r' * r) / (oldr' * oldr);
        d = r + beta * d;
    end
end

