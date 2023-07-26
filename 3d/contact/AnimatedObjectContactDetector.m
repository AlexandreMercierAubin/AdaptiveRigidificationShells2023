classdef AnimatedObjectContactDetector < ContactFinder
    %ObjectContactDetector collision detector that allows the collision
    %detection between meshes and an imported obj mesh
    
    properties
        parts
        skipdetection = false;
        interpolateOnFrames = 1;
        plotScale = 1;
        stopFrame = 0;
        edgeColor = 'none';
        F
        faceLists
        h
        sizeV= [0,0];
        animationPositionStart = 0;
        skipPart = [];
        holdPose = [];%maybe make it logical for quick lookup
        frameHeld = 0;
    end
    
    methods
        function obj = AnimatedObjectContactDetector(filename, partsFolder, degrees, scale, position, frictionCoefficient, h)
            disp("loading animated collider");
            [positions,orientations,parts] = readMOT(filename);
            
            obj.parts = AnimatedPart.empty;

            angles = deg2rad(degrees);
            quatAngles = eul2quat(angles);
            R = quat2mat(quatAngles);
            
            for i = 1:numel(parts)
                scaledPositions = scale.*positions(:,:,i);
                transformedPositions = (R*scaledPositions')';
                %matlab quaternions have the first and last entries inverted
                rotatedQuats = quatmultiply(quatAngles,orientations(:,:,i));
                rotatedQuats = normr(rotatedQuats);
                obj.parts = [obj.parts, AnimatedPart(partsFolder,parts{i},transformedPositions , rotatedQuats, scale, position)];
                obj.sizeV = [obj.sizeV(1)+size(obj.parts(i).VerticesRest,1),3];
            end
            obj.skipPart = false(numel(parts),1);
            
            obj.FrictionCoefficient = frictionCoefficient;

            Fcells = {obj.parts.Faces};
            obj.faceLists = zeros(numel(obj.parts)+1,1);
            obj.F = Fcells{1};
            for i = 2:numel(obj.parts)
                maxF = max(obj.F(:));
                obj.faceLists(i) = maxF;
                obj.F = [obj.F;Fcells{i}+maxF];
            end
            obj.faceLists(numel(obj.faceLists)) = max(obj.F(:));
            obj.h = h;
            disp("done loading animated collider");
        end
        
        function obj = injectPause(obj,timeArray,pauseDuration)
            totalFrames = obj.interpolateOnFrames * size(obj.parts(1).Orientations,1);
            obj.holdPose = false(totalFrames,1);
            for i = numel(timeArray):-1:1 %last to first so the indexes remain valid
                pos = timeArray(i);
                trueVec = true(pauseDuration(i),1);
                obj.holdPose = [obj.holdPose(1:pos);trueVec;obj.holdPose((pos+1):end)];
            end
        end

        function [Jc, phi, cInfo] = findContacts( obj, meshes, time, p)
            phi = zeros( 0, 1 );
            Jc = zeros( 0, meshes.N*3 );
            cInfo = contactInfo3D.empty;
            obj.interpolateOnFrames = max(obj.interpolateOnFrames,1);
            
            animFrame = time - obj.frameHeld;
            useTrajectory = true;
            if ~isempty(obj.holdPose) && time <= numel(obj.holdPose) && obj.holdPose(time)
                useTrajectory = false;
                obj.frameHeld = obj.frameHeld + 1;
            end

            for i = 1:numel(obj.parts)
                if obj.skipPart(i)
                    continue;
                end
                [JcOut, phiOut, cInfoOut] = obj.findPartContact(meshes,animFrame,p,obj.parts(i),useTrajectory);
                if ( isempty(phiOut) )
                    continue;
                end
                Jc = [ Jc; JcOut ];
                phi = [ phi; phiOut ];
                cInfo = [ cInfo, cInfoOut ];
            end
        end
        
        function [J, phi, cInfo]= findPartContact(obj,meshes,time,p,part, useTrajectory)
            ps = vertcat(p);
            xs = ps(1:3:end);
            ys = ps(2:3:end);
            zs = ps(3:3:end);
            
            [V] = obj.getObjPositionPart(time,1,part);
            
            % AABB
            [minV,maxV] = bounds(V,1);
            AABB = [minV',maxV'];
            vertThickness = meshes.vertexContactThickness;
            s1InAABB2 = ...
                xs + vertThickness >= AABB(1,1) & ...
                xs - vertThickness <= AABB(1,2) & ...
                ys + vertThickness >= AABB(2,1) & ...
                ys - vertThickness <= AABB(2,2) & ...
                zs + vertThickness >= AABB(3,1) & ...
                zs - vertThickness <= AABB(3,2) ;
            idx = s1InAABB2;

            if ( sum(idx) == 0  || obj.skipdetection)
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            points = [xs(idx),ys(idx),zs(idx)];
            
           [S,I,C,N] = signed_distance(points, V, part.Faces);
        
            % these interpenetration depths must be negative for
            % baumgarte to work!
            indices = find(idx);
            hitThickness = meshes.vertexContactThickness(indices);
            
            in = S < hitThickness & ~isnan(N(:,1)) & ~isnan(N(:,2)) & ~isnan(N(:,3));
            
            points = points(in,:);
            
            if sum(in) <= 0
                J = zeros(0,numel(ps));
                phi = [];
                cInfo = contactInfo3D.empty;
                return;
            end
            
            phi = S(in) - hitThickness(in);
            phi(phi > 0) = 0;
            
            if useTrajectory
                CV = obj.getTrajectory(time, points, part);
            else
                CV = zeros(size(points));
            end

            normal = N(in,:);
            tangent = zeros(size(normal));
            tangent2 = zeros(size(normal));
            for i = 1: size(normal,1)
                tangents = null(normal(i,:));
                tangent(i,:) = tangents(:,1)';
                tangent2(i,:) = tangents(:,2)';
            end
            
            indices = indices(in);

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
            for i = 1:numel(phi) 
                cInfo(i) = contactInfo3D( [xs(indices(i)), ys(indices(i)), zs(indices(i))], normal(i,:), tangent(i,:), obj.FrictionCoefficient, indices(i), obj.ID, tangent2(i,:));
                cInfo(i).velocity = -[dot(CV(i,:)',normal(i,:)),dot(CV(i,:)',tangent(i,:)),dot(CV(i,:)',tangent2(i,:))];
            end             
        end

        function render( obj, frame )
            [V] = obj.getObjPositionFaces(frame, obj.plotScale);
            
            if ( obj.plotHandle ~= 0 )
                obj.plotHandle.Vertices = V;
                return
            end

            hold on;
            obj.plotHandle = patch('Faces',obj.F,'Vertices',V,'FaceColor','red', 'EdgeColor', obj.edgeColor);
        end
        
        function [CV] = getTrajectory(obj,frame, P, part)
            CV = zeros(size(P));
            [time, alpha] = computeTimeAlpha(obj,frame);
            [time2, alpha2] = computeTimeAlpha(obj,frame+1);
            if (obj.stopFrame ~= 0 && frame > obj.stopFrame) || obj.parts(1).isAnimationOver(time)
                return;
            end

            points = [P,ones(size(P,1),1)];

            M0inv = part.getInverseTransform(time,1);
            M1inv = part.getInverseTransform(time+1,1);
            Minv = alpha*M1inv + (1-alpha)*M0inv;

            M0 = part.getTransform(time2,1);
            M1 = part.getTransform(time2+1,1);
            M = alpha2*M1 + (1-alpha2)*M0;

            Mchange = M*Minv;
            nextPoints = (Mchange*points')';
            nextPoints = nextPoints(:,1:3)./nextPoints(:,4);
            CV =  (nextPoints-P)./obj.h;
        end

        function [time, alpha] = computeTimeAlpha(obj,frame)
            time = floor(frame/obj.interpolateOnFrames)+1 + obj.animationPositionStart; %TODO: make sure the number of divisions fits otherwise it might cause awkward pauses
            modValue = mod(frame,obj.interpolateOnFrames);
            alpha = modValue / obj.interpolateOnFrames;
        end


        function [CV] = getTrajectories(obj,frame,P,I)

            CV = zeros(size(P));
            [time, alpha] = computeTimeAlpha(obj,frame);
            [time2, alpha2] = computeTimeAlpha(obj,frame+1);
            if (obj.stopFrame ~= 0 && frame > obj.stopFrame) || obj.parts(1).isAnimationOver(time)
                return;
            end

            for i = 2:numel(obj.parts)+1
                partID = i-1;
                pointLogical = I<= obj.faceLists(i) & I > obj.faceLists(i-1);
                if(~any(pointLogical))
                    continue;
                end
                points = [P(pointLogical,:),ones(sum(pointLogical),1)];

                M0inv = obj.parts(partID).getInverseTransform(time,1);
                M1inv = obj.parts(partID).getInverseTransform(time+1,1);
                Minv = alpha*M1inv + (1-alpha)*M0inv;

                M0 = obj.parts(partID).getTransform(time2,1);
                M1 = obj.parts(partID).getTransform(time2+1,1);
                M = alpha2*M1 + (1-alpha2)*M0;

                Mchange = M*Minv;
                nextPoints = (Mchange*points')';
                nextPoints = nextPoints(:,1:3)./nextPoints(:,4);
                CV(pointLogical,:) =  (nextPoints-points(:,1:3))./obj.h;
            end

        end

        function [V] = getObjPositionPart(obj, frame, scaleRatio, part)
            if obj.stopFrame ~= 0 && frame > obj.stopFrame
                frame = obj.stopFrame;
            end
            
            [time, alpha] = computeTimeAlpha(obj,frame);

            if frame == 1
                M = part.getTransform(time,scaleRatio);
            else
                M0 = part.getTransform(time,scaleRatio);
                M1 = part.getTransform(time+1,scaleRatio);
                
                M = alpha*M1 + (1-alpha)*M0;
            end
            V = (M * [part.VerticesRest,ones(size(part.VerticesRest,1),1)]')';
            V = V(:,1:3) ./ V(:,4);
        end

        function [V,F] = getObjPositionFaces(obj, frame, scaleRatio)
            if nargin < 3
                scaleRatio = obj.plotScale;
            end
            animFrame = frame - obj.frameHeld;

            V = zeros(obj.sizeV); %TODO: there should always be the same number of verts so make this a constant size
            if obj.stopFrame ~= 0 && animFrame > obj.stopFrame
                animFrame = obj.stopFrame;
            end
            
            [time, alpha] = computeTimeAlpha(obj,animFrame);
            counter = 1;

            for i = 1:numel(obj.parts)
                if animFrame == 1
                    M = obj.parts(i).getTransform(time,scaleRatio);
                else
                    M0 = obj.parts(i).getTransform(time,scaleRatio);
                    M1 = obj.parts(i).getTransform(time+1,scaleRatio);
                    
                    M = alpha*M1 + (1-alpha)*M0;
                end
                Vobject = (M * [obj.parts(i).VerticesRest,ones(size(obj.parts(i).VerticesRest,1),1)]')';
                Vobject = Vobject(:,1:3) ./ Vobject(:,4);
                V(counter:counter+size(Vobject,1)-1,:) = Vobject;
                counter = counter + size(Vobject,1);
            end
            if nargout >1
                F = obj.F;
            end
        end
    end
end

