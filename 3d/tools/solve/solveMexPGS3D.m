function [lambda] = solveMexPGS3D(iterations, A, b, initial_lambda, cInfo, compliance, td )
    % SOLVECLDLPGS PGS solver for 2D contact constraints that uses the LDLT
    % factorization of A = M - h^2K, but does a full assembly of the JcAinvJT matrix 

    n = size(b, 1);
    lambda = initial_lambda;
    if ~n
        return;
    end 
            
    mu = [ cInfo(:).frictionCoefficient ];
    [ lambda ] = mexPGS3DwithJAinvJT( iterations, lambda, A, b, mu, compliance );
%     [ lambda ] = solveLDLTPGS3DwithJAinvJTHelper( iterations, lambda, A, b, mu, compliance );
    
end
