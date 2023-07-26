function [F] = computeFtri(V,T,DmInv)
    numt = size(T,1);
    F = zeros(numt*4,1);
    for i = 1:numt
            %Fetching the position of eatch vertex
            t = T(i,:);
            p0 = V(t(1),:);
            p1 = V(t(2),:);
            p2 = V(t(3),:);
            
            %femdefo's Dm matrix, this is essentially a material frame of the
            %undeformed mesh
            localT = [p1-p0;p2-p0];
            Ds = localT*localT';
            F(i*4-3:i*4) = Ds*DmInv(:,:,i); 
    end
end

