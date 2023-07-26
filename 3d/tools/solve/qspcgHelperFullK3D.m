function Av = qspcgHelperFullK3D( h, M, D, K, v )
% QSPCGHELPER Provides a multiplication by system matrix A without assembly
% and does so with an elimination of common subexpressions
    Av = M*v - h*(D*v) - (h*h)*(K*v);
end