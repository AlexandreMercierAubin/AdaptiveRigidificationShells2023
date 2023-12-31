function AinvBlocks = computeAInverseDiagonalBlocks3D(A)
%COMPUTEAINVERSEDIAGONALBLOCKS Computes a sparse matrix with 3x3 blocks of
%only the diagonal of A inverse.
%   There seems to be no really good efficient way to do this.  Probably a
%   good idea to cache the answer.  Could be make slightly faster with a
%   selective forward substitution.

    n = size(A,1);
    [L, D, P, S] = ldl(A);
    iD = 1 ./ full(diag(D)); % not sure this is faster than backslash??  was this tested?
    
    numBlocks = n/3;    % n must be invisible by 3
    blocks = cell(numBlocks,1);
    parfor block = 1:numBlocks 
        binds = (block-1)*3+1:block*3;
        b = zeros( n, 3 );
        b(binds,:) = eye(3);
        x = S * (P * (L' \ ( iD .* (L \ (P' * (S * b))))));
        blocks{block} = x(binds,:);
    end
    AinvBlocks = blkdiag( sparse(blocks{1}), blocks{2:end} );
end

