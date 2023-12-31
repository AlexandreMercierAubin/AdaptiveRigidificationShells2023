classdef EnergyModel < handle
    %EnergyModel abstract class that allows energies to be computed modularly
    properties
        energy
        derivative1Gradient
        derivative2HessianC
        elasticForces
        stiffness
        stiffnessDamping
        name = "default"
    end
    
    methods
        function obj = EnergyModel()
            obj.energy = NaN; 
            obj.derivative1Gradient = NaN;
            obj.derivative2HessianC = NaN;
            obj.elasticForces = NaN;
            obj.stiffness = NaN;
            obj.stiffnessDamping = NaN;
        end 
    end
    methods (Abstract = true)
            computeEnergy(obj, mesh3D, deformationGradientF);
    end
end

