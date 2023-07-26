function plotDebugSpheresAtPos(pos)
% plotDebugSpheresAtPos(pos)
% pos : n by 3 matrix of positions x,y,z
    for i = 1 : size(pos,1)
        c = pos(i,:);
        r = 0.01;
    
        [x,y,z] = sphere;
        x = x * r;
        y = y * r;
        z = z * r;
    
        xdata = c(1)+x;
        ydata = c(2)+y;
        zdata = c(3)+z;
        
        cdata = ones(size(zdata));
        cdata(1) = 0;
        
        hold on;
        plotHandle = surf(xdata,ydata,zdata,cdata,'EdgeColor','none','FaceColor',[0 0 0]);
    end
end

