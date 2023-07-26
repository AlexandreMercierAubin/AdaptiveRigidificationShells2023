function x = solveCGNoLoop(A, b, initial_x)
    x = initial_x;
    r = b - A(x);
    alpha = (r' * r) / (r' * A(r));
    x = x + alpha * r;
end

