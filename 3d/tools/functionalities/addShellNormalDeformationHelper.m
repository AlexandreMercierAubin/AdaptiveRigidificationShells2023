function [n_vec, ii,jj,vals, valCount] = addShellNormalDeformationHelper(V, T, isShellElement, referenceSpaceNormals, normals)
        % we need to add the normal change to the deformation gradient of shells
        numt = size(T,1);
        n_vec = zeros(9*numt,1);
%         cache.dndp = zeros(3,9,numt);
        localN = zeros(9,3);

        %preallocates the sparse matrix with the biggest case scenario in
        %mind (all shells)... can be improved with the number of shell
        %elements instead
        ii = zeros(1,9*9*numt);
        jj = zeros(1,9*9*numt);
        vals = zeros(1,9*9*numt);
        counter = 0;

        for i = 1:numt
            t = i;
            if ~isShellElement(t)
                continue;
            end

            v1Id = T(t,1);
            p0 = V(v1Id,:)';
            v2Id = T(t,2);
            p1 = V(v2Id,:)';
            v3Id = T(t,3);
            p2 = V(v3Id,:)';
            
            localN(1:3,1) = referenceSpaceNormals(i,:);
            localN(4:6,2) = referenceSpaceNormals(i,:);
            localN(7:9,3) = referenceSpaceNormals(i,:);

            localnvec =  localN*normals(i,:)';
            n_vec(t*9-8:t*9) = localnvec;

            if nargout == 1
                continue;
            end

            x = [p0;p1;p2];
            I1 = [-eye(3), eye(3), zeros(3)];
            I2 = [-eye(3), zeros(3), eye(3)];
            laplacex1 = I1 * x;
            laplacex2 = I2 * x;

            ntilde = cross(laplacex1,laplacex2);

            crossx1 = crossProductMatrix(laplacex1);
            crossx2 = crossProductMatrix(laplacex2);
            
            normntilde = norm(ntilde);
            assert(~isinf(normntilde));%if this triggers, then the normal is made of zeros... check your mesh
            smalln = ntilde./normntilde;          
            
            invNtilde = 1/normntilde;

            dndp = invNtilde*(eye(3)-smalln*smalln')*(crossx1*I2 - crossx2*I1); 
%             cache.dndp(:,:,t) = dndp;
            Ndndp = localN*dndp;

            counter = counter + 1;
            idi = repmat(t*9-8:t*9,1,9);
            ii(counter*81-80:counter*81) = idi;
            idj1 = reshape(repmat(T(t,1)*3-2:T(t,1)*3,9,1),[],1)';
            idj2 = reshape(repmat(T(t,2)*3-2:T(t,2)*3,9,1),[],1)';
            idj3 = reshape(repmat(T(t,3)*3-2:T(t,3)*3,9,1),[],1)';
            jj(counter*81-80:counter*81) = [idj1,idj2,idj3];
            B1 = reshape(Ndndp(1:9,1:3),1,[]);
            B2 = reshape(Ndndp(1:9,4:6),1,[]);
            B3 = reshape(Ndndp(1:9,7:9),1,[]);
            vals(counter*81-80:counter*81) = [B1,B2,B3];
        end
        valCount = counter*81;
end