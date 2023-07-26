function Edot = strainRate(F,Fprev,h,dim)
    s = strain(F,dim);
    sprev = strain(Fprev,dim);
    Edot = (s-sprev)./h;
end

