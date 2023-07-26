function [E] = strain(F,dim)
    if dim ==3
        n = 9;
    else
        n = 4;
    end
    E = zeros(size(F));
    for i = 1:numel(F)/n
        ids = (i*n-(n-1)):i*n;
        Fvec = F(ids);
        Flocal = reshape(Fvec,dim,dim);
        E(ids) = 0.5*(Flocal'*Flocal - eye(dim));
    end
end

