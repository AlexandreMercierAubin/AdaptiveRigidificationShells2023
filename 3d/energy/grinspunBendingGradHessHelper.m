function [ii,jj, vals, valsD, fullg]  = grinspunBendingGradHessHelper(p, edgeNormalConsistentVertices, edgeOppositeVertices, kd, anglesDiff, norm1, norm2, area, areaTilde, innerAngle1, innerAngleTilde1, innerAngle2, innerAngleTilde2, bendingEdgeAlpha1)
        fullg = zeros(size(p,1),1);
        numEdges = size(edgeNormalConsistentVertices,1);
        ii = zeros(numEdges*144,1);
        jj = zeros(numEdges*144,1);
        vals = zeros(numEdges*144,1);
        valsD = zeros(numEdges*144,1);
        for i = 1:numEdges
            vertex1 = edgeNormalConsistentVertices(i,1);
            vertex2 = edgeNormalConsistentVertices(i,2);
            opposite1 = edgeOppositeVertices(i,1);
            opposite2 = edgeOppositeVertices(i,2);

            v1 = p(vertex1*3-2:vertex1*3);
            v2 = p(vertex2*3-2:vertex2*3);
            o1 = p(opposite1*3-2:opposite1*3);
            o2 = p(opposite2*3-2:opposite2*3);

            ids = [opposite1*3-2:opposite1*3,vertex1*3-2:vertex1*3,vertex2*3-2:vertex2*3,opposite2*3-2:opposite2*3];
            
            [g,H] = neoGrinspunianBendingEnergy(kd(i),anglesDiff(i),[o1';v1';v2';o2'],norm1(i,:),norm2(i,:),area(i),areaTilde(i),innerAngle1(i),innerAngleTilde1(i),innerAngle2(i),innerAngleTilde2(i));

            idi = repmat(ids',1,12);
            ii(i*144-143:i*144) = idi;
            idj1 = reshape(repmat(ids,12,1),[],1)';
            jj(i*144-143:i*144) = idj1;
            vals(i*144-143:i*144) = H(:);
            valsD(i*144-143:i*144) = bendingEdgeAlpha1(i)*H(:);

            fullg(ids) = fullg(ids) + g;

        end
end

