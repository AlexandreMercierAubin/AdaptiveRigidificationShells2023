classdef RigidBody3D < handle
    properties
        Mesh            % Mesh of this RigidBody
        
        Position        % rigid position
        Velocity        % rigid body velocity
        Rotation        % rigid orientation
        AngularVelocity % rigid angular velocity

        Force           % rigid force accumulator
        Torque          % rigid torque accumulator

        Mass            % mass of rigid
        Inertia         % inertia of the rigid body
        Inertia0        % body inertia at identity rotation

        DOFs            % list of rigid DOFs in the 2 * N state vector
        VertexIDs         % indices of rigid vertices 
        TetInds         % rigid tetrahedron indices in this body
        isPinned
        isAnimated
    end
    
    methods
        function obj = RigidBody3D(mesh)
            %RIGIDBODY makes a rigid body to be part of an AdaptiveMesh.
            %   Empty rigid bodies should not exist (rigid bodies with no
            %   tet)
            obj.Mesh = mesh;
            obj.Position = zeros(3, 1);
            obj.Velocity = zeros(3, 1);
            obj.Rotation = eye(3);

            obj.Force = zeros(3, 1);
            obj.Torque = zeros(3, 1);

            obj.Mass = 0;
            obj.Inertia = zeros(3,3); 
            obj.Inertia0 = zeros(3,3); 
            
            obj.DOFs = [];
            obj.VertexIDs = [];
            obj.TetInds = [];
            obj.AngularVelocity = zeros(3,1);
            obj.isPinned = false;
        end
        
        function clone = clone(obj, mesh)
            clone = RigidBody3D();
            fns = properties(obj);
            for i = 1:numel(fns)
                clone.(fns{i}) = obj.(fns{i});
            end
            clone.Mesh = mesh;
        end
        
        function addVertices(obj, indices, updateMesh)
            if nargin < 3
                updateMesh = 1;
            end
            
            obj.VertexIDs = unique([obj.VertexIDs, indices]);
            obj.updateBody(updateMesh);
        end
        
        function removeVertices(obj, indices)
            obj.VertexIDs(ismember(indices, obj.VertexIDs)) = [];
            obj.updateBody();
        end
        
        removeTriangles(obj, tris, updateMesh, settings);
        
        function merge(obj, body, updateMesh)
            if nargin < 3
                updateMesh = 1;
            end
            obj.Mesh.RigidBodies(obj.Mesh.RigidBodies == body) = [];
            obj.VertexIDs = unique([obj.VertexIDs, body.VertexIDs]);
            obj.updateBody(updateMesh);
        end
        
        function updatePosition(body, h, deltav)
            % updatePosition Performs the symplectic update of velocity and
            % the position for both linear and angular velocity
            body.Velocity = body.Velocity + deltav(1:3);
            body.Position = body.Position + h * body.Velocity;
            if body.isAnimated
                %prevent rotation for scripted animations
                body.AngularVelocity = body.AngularVelocity;
            else
                body.AngularVelocity = body.AngularVelocity + deltav(4:6);
            end
            body.Rotation = expRodrigues( body.AngularVelocity, h ) * body.Rotation;
            body.Inertia = body.Rotation * body.Inertia0 * body.Rotation'; 
            body.setDOFsFromRigid(h);
        end

        function [pList, vList,inds, bodyp, angularv,inertia] = peekPosition(body, h, deltav, mesh3Dp)
            % updatePosition Performs the symplectic update of velocity and
            % the position for both linear and angular velocity
            bodyp = body.Position + h * deltav(1:3);
            if body.isAnimated
                %prevent rotation for scripted animations
                angularv = body.AngularVelocity;
            else
                angularv = body.AngularVelocity + deltav(4:6);
            end
            R = expRodrigues( angularv, h ) * body.Rotation;
            inertia = body.Rotation * body.Inertia0 * body.Rotation';
            [pList, vList, inds] = body.peekDOFsFromRigid(h, R, bodyp, mesh3Dp);
        end

       
        %TODO: move this to the rigidbody update thingy so it is only computed
        %when rigid bodies are changed. This function will only extract the
        %values of alpha0
        function alpha0 = getAlpha0(body, gamma)
            ids = body.VertexIDs * 3;
            alpha0s = zeros(size(ids,2)*3,1);
            alpha0s(3:3:end) = body.Mesh.alpha0(ids);
            alpha0s(2:3:end) = body.Mesh.alpha0(ids-1);
            alpha0s(1:3:end) = body.Mesh.alpha0(ids-2);
            croppedGamma = zeros(size(ids,2)*3,3);
            croppedGamma(3:3:end) = gamma(ids,:);
            croppedGamma(2:3:end) = gamma(ids-1,:);
            croppedGamma(1:3:end) = gamma(ids-2,:);
            alpha0 = zeros(3,1);
            
            for j = 1:1:size(gamma,2)
                alpha0(j) = sum(alpha0s.*croppedGamma(:,j));
            end
        end
        
        computeRigidForce( body, h);
        [Force, Torque] = peekRigidForce(body, h, AngularVelocity, RigidBodyPosition, Inertia, meshf, meshp);
        updateBody( obj, updateMesh );  
        setDOFsFromRigid( body, h );
    end
end

