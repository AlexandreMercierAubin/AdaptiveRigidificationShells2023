function x = solvePCG( iterations, A, b, initial_x, cache, resStop)
% SOLVEPCG Conjugate Gradient on Ax=b, using a preconditioner and starting 
% with the initial point (should be zero! Use shifted Krylov to warm start
% instead).  
%   iterations  maximum number of iteraitons
%   A           matlab function that multiplies by A
%   b           the right hand side
%   initial_x   should be zero!  TODO: remove?
%   cache       has the preconditioner, a matlab function that applies the 
%               preconditioner. TODO: pass only the preconditioner?    x = initial_x;
    if nargin < 6
        resStop = 1e-30;
    end
    x = initial_x;
    r = b - A(x);
    z = cache.preconditioner(r);
    q = z;
    rTz = r' * z;
    for i=1:iterations
        Ap = A(q);
        qAp = (q'*Ap);
        if qAp == 0
            alpha = 0;
        else
            alpha = rTz / (q'*Ap);
        end
        x = x + alpha * q;
        r = r - alpha * Ap;
        if norm(r) < resStop
              break;
        end
        z = cache.preconditioner(r);
        rTznew = r'*z;
        beta = rTznew / rTz;
        rTz = rTznew;
        q = z + beta * q;
    end
end