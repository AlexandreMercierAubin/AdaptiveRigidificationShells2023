function [plotSkip] = plotSkip60FPS(h)
%PLOTSKIP60FPS compute the plotskip needed for 60 rendering per second
%according to substep h
FPS = 1/60;    
scale = 1/h *FPS;
plotSkip = floor(scale);
end

