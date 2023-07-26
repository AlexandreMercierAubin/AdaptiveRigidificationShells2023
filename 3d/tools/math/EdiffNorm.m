function [Edn] = EdiffNorm(F,Fprev,h,dim)
    if dim ==3
        n = 9;
    else
        n = 4;
    end
    Edot = strainRate(F,Fprev,h,dim);
    EdnSize = numel(F)/n;
    Edn = zeros(EdnSize,1);
    for i = 1 : EdnSize
        ids = (i*n-(n-1)):i*n;
        EdotLocal = reshape(Edot(ids),dim,dim);
        Edn(i) = norm(EdotLocal,"fro")^2;
    end
end

