function res = mult2x2Stack(A, B)
    a = A(1, 1, :) .* B(1, 1, :) + A(1, 2, :) .* B(2, 1, :);
    b = A(1, 1, :) .* B(1, 2, :) + A(1, 2, :) .* B(2, 2, :);
    c = A(2, 1, :) .* B(1, 1, :) + A(2, 2, :) .* B(2, 1, :);
    d = A(2, 1, :) .* B(1, 2, :) + A(2, 2, :) .* B(2, 2, :);

    res(1:2, 1:2, :) = [a, b; c, d];
end