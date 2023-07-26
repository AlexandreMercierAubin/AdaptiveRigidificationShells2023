classdef TriangleMaterial < handle
    %Material that can independently set for each element, and/or reused for
    %multiple of them.
    properties
        rho = 10;     % Density
        mu      % Lamé parameter
        lambda  % Lamé parameter
        alpha0 = 0.03; % Rayleigh parameter
        alpha1 = 0.02;  % Rayleigh parameter
        bendingAlpha1 = 1;
        strainUpperBound = Inf; % 5% would be 1.05
        strainLowerBound = -Inf;% 5% would be 0.95
        color = [0.3,0.46,0.8]           % color of the material
        cacheName = 'defaultCacheName';
        thickness = 0.1;
        grinspunShearingStiffness = 1;
        grinspunAreaStiffness = 0.5;
    end
    
    methods
        function obj = TriangleMaterial(rho, mu, lambda, alpha0, alpha1, color, strainUpperBound, strainLowerBound, thickness, kl, ka, bendingAlpha1)
            if nargin > 0
                obj.rho = rho;
            end
            
            if nargin > 1
                obj.mu = mu;
            else
                k = 4e4;
                nu = 0.40;
                mu = k / 2 / (1 - nu); 
                obj.mu = mu;
            end
            
            if nargin > 2
                obj.lambda = lambda;
            else
                k = 4e4;
                nu = 0.40;
                lambda = k * nu / (1 + nu) / (1 - 2 * nu);
                obj.lambda = lambda;
            end
            
            if nargin > 3
                obj.alpha0 = alpha0;
            end
            
            if nargin > 4
                obj.alpha1 = alpha1;
            end

            if nargin > 5
                obj.color = color;
            end
            
            if nargin > 6
                obj.strainUpperBound = strainUpperBound;
            end
            
            if nargin > 7
                obj.strainLowerBound = strainLowerBound;
            end

            if nargin > 8
                obj.thickness = thickness;
            end

            if nargin > 9
                obj.grinspunShearingStiffness = kl;
            end

            if nargin > 10
                obj.grinspunAreaStiffness = ka;
            end

            if nargin > 11
                obj.bendingAlpha1 = bendingAlpha1;
            else
                obj.bendingAlpha1 = obj.alpha1;
            end
        end
    end
end
