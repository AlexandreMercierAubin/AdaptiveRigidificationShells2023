classdef MovingSphereContactFinder < ContactFinder
    %WIP SphereContactFinder Contact finder class that only computes contact
    %between meshes and a sphere defined by normal and position. The sphere
    %moves according to a callback function
    
    properties
        Radius
        Center
        
        % NOTE: trajectory could be a property rather than hard coded...
        % could be a function that could be evaluated.
        cfun
        dcdt
        thetafun
        dthetadt
        plotRatio
        V
        F
        h
    end
    
    methods
        function obj = MovingSphereContactFinder( radius, position, frictionCoefficient,h )
            obj.h = h;
            obj.Radius = radius;
            obj.Center = position;
            if nargin >= 3
                obj.FrictionCoefficient = frictionCoefficient;
            end
            
            obj.cfun = @(t) [ 0, 0.1*(cos(t*8)-1) * (mod( t*8, pi*2*10 ) < pi*2), 0];
            obj.dcdt = @(t) [ 0, -0.1*sin(t*8)*8  * (mod( t*8, pi*2*10 ) < pi*2), 0];

            obj.thetafun = @(t)  0.03*sin(t*6)    * (mod( t*6, pi*2*3 ) < pi*2);
            obj.dthetadt = @(t)  0.03*cos(t*6)*4  * (mod( t*6, pi*2*3 ) < pi*2);
            obj.plotRatio = 1.0;
            [obj.V,obj.F] = subdivided_sphere(3);
        end
        
        function [J, phi, cInfo] = findContacts( obj, meshes, time , p)
            ps = vertcat(p);
            V = meshes.formatPositions(p);
            
            c = obj.Center + obj.cfun(time);
            distance = vecnorm(c-V,2,2);
            
            % indices of the ones under the plane
            idx = distance <= obj.Radius;
            
            if ( all(idx == 0) )
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            % constraint value (penetration)
            phi = distance(idx) - obj.Radius;
            normal = [V(idx,1) - c(1), V(idx,2) - c(2), V(idx,3) - c(3)] ./ distance(idx);
            tangent = zeros(size(normal));
            tangent2 = zeros(size(normal));
            for i = 1: size(normal,1)
                tangents = null(normal(i,:));
                tangent(i,:) = tangents(:,1)';
                tangent2(i,:) = tangents(:,2)';
            end
            
            indices = find(idx);

            %add friction
            rown  = (1:3:3*numel(phi))';
            rowt  = (2:3:3*numel(phi))';
            rowt2  = (3:3:3*numel(phi))';
            colx = (indices*3-2);
            coly = (indices*3-1);
            colz = (indices*3);
            J = sparse( ...
                [ rown;         rown;      rown;       rowt;         rowt;     rowt;   rowt2;         rowt2;     rowt2], ...
                [ colx;         coly;      colz;       colx;         coly;     colz;   colx;         coly;     colz;], ...
                [ normal(:,1);  normal(:,2); normal(:,3); tangent(:,1); tangent(:,2) ; tangent(:,3) ; tangent2(:,1); tangent2(:,2) ; tangent2(:,3)], 3*numel(phi), numel(ps) );

            cInfo = contactInfo3D.empty; %cell( numel(phi), 1 );
            R = eul2rotm([deg2rad(obj.h*45),0,0]);%obj.dthetadt(time)
            rotatedContacts = (R*(V(indices,:)-c)')'+c;
            vel = (rotatedContacts - V(indices,:))./obj.h;
            for i = 1:numel(phi) 
                cInfo(i) = contactInfo3D( [V(indices(i),1), V(indices(i),2), V(indices(i),3)], normal(i,:), tangent(i,:), obj.FrictionCoefficient, indices(i), obj.ID, tangent2(i,:));
                cInfo(i).velocity = -[dot(normal(i,:),obj.dcdt(time)); dot(tangent(i,:),obj.dcdt(time))- dot(tangent(i,:),vel(i,:)); dot(tangent2(i,:),obj.dcdt(time))- dot(tangent(i,:),vel(i,:))]; 
            end             
        end
        
        function render( obj, time)
            if ( obj.plotHandle ~= 0 )
                return
            end
            hold on;
            V = obj.Center + obj.cfun(time) + obj.V*obj.Radius*obj.plotRatio;
            obj.plotHandle = patch('Faces',obj.F,'Vertices',V,'FaceColor','Blue','EdgeColor','none');
        end
        
        function [V,F] = getObjPositionFaces(obj, time)
            F = obj.F;
            V = obj.Center + obj.V*obj.Radius*obj.plotRatio;
        end
    end
end

