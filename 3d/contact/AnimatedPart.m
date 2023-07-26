classdef AnimatedPart < handle
    %ANIMATEDPART Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        VerticesRest
        Faces
        Positions    %n by 3 position of the reference frame
        Orientations %n by 4 quaternion orientation of the part
        Offset % starting position
    end
    
    methods
        function obj = AnimatedPart(infolder, inentry, positions, orientations, scale, offset)
            filename = [infolder,inentry{1},num2str(inentry{2},'%04.f'),'.obj'];
            [obj.VerticesRest,obj.Faces,UV,TF,N,NF] = readOBJ(filename);
            obj.VerticesRest = scale.*obj.VerticesRest;
            obj.Positions = positions;
            obj.Orientations = orientations;
            obj.Offset = offset;
        end

        function isOver = isAnimationOver(obj,frame)
            isOver = frame > size(obj.Orientations,1);
        end

        function M = getInverseTransform(obj, frame, scaleRatio)
            refFrame = frame;
            if isAnimationOver(obj,frame)
                refFrame = size(obj.Orientations,1);
            end
            s = 1/scaleRatio;
            R = eye(4);
            R(1:3,1:3) = quat2mat(obj.Orientations(refFrame,:));
            S = diag([s,s,s,1]);
            T = eye(4);
            T(1:3,4) = -(obj.Positions(refFrame,:) + obj.Offset);
            M = S*R'*T;
        end

        function vertices = getVerticesAtFrame(obj, frame, scaleRatio)
            refFrame = frame;
            if isAnimationOver(obj,frame)
                refFrame = size(obj.Orientations,1);
            end
            
            R = quat2mat(obj.Orientations(refFrame,:));

            rotatedVerts = (R*(scaleRatio.*obj.VerticesRest)')';
            vertices = rotatedVerts + obj.Positions(refFrame,:) + obj.Offset;

        end

        function M = getTransform(obj, frame, scaleRatio)
            refFrame = frame;
            if frame > size(obj.Orientations,1)
                refFrame = size(obj.Orientations,1);
            end
            R = eye(4);
            R(1:3,1:3) = quat2mat(obj.Orientations(refFrame,:));
            S = diag([scaleRatio,scaleRatio,scaleRatio,1]);
            T = eye(4);
            T(1:3,4) = obj.Positions(refFrame,:) + obj.Offset;
            M = T*R*S;
        end
        
    end
end

