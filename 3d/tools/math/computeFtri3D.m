function [F] = computeFtri3D(V,T,N,dphidx)
    assert(size(T,1) == size(N,1));
    numt = size(T,1);
    F = zeros(numt*4,1);
    for i = 1:numt
            %Fetching the position of eatch vertex
            t = T(i,:);
            p0 = V(t(1),:);
            p1 = V(t(2),:);
            p2 = V(t(3),:);
            n = N(i,:);
            
            F(i*9-8:i*9) = [p0;p1;p2;n]'*dphidx(:,:,i); 
    end
end

