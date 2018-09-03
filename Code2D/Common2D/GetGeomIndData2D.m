fprintf('... generating geometric data needed by indicators\n')

% Find neighbors in patch
E1 = EToE(:,1)'; E2 = EToE(:,2)'; E3 = EToE(:,3)';

% Find local faces in neighbouring triangles
f1 = EToF(:,1)'; f2 = EToF(:,2)'; f3 = EToF(:,3)';

% Extract coordinates of vertices and centers of elements
v1 = EToV(:,1); xv1 = VX(v1); yv1 = VY(v1);
v2 = EToV(:,2); xv2 = VX(v2); yv2 = VY(v2);
v3 = EToV(:,3); xv3 = VX(v3); yv3 = VY(v3);  


% Get vertices opposite shared face in neighbour
vind  = [3,1,2];
xvn1  = zeros(1,K); xvn2 = xvn1; xvn3 = xvn1;
yvn1  = zeros(1,K); yvn2 = yvn1; yvn3 = yvn1;
for i = 1:K
   vn1   = EToV(E1(1,i),vind(f1(1,i))); xvn1(1,i) = VX(vn1); yvn1(1,i) = VY(vn1);
   vn2   = EToV(E2(1,i),vind(f2(1,i))); xvn2(1,i) = VX(vn2); yvn2(1,i) = VY(vn2);
   vn3   = EToV(E3(1,i),vind(f3(1,i))); xvn3(1,i) = VX(vn3); yvn3(1,i) = VY(vn3);
   
   % Correctional shift periodic ghost triangles
   if(isKey(PShift,i) && UseMeshPerData)
        pairdata = PShift(i);
        p_ind    = find(pairdata(:,1)==E1(1,i));
        if ~isempty(p_ind)
            xvn1(1,i) = xvn1(1,i) - pairdata(p_ind,2);
            yvn1(1,i) = yvn1(1,i) - pairdata(p_ind,3);
        end
        p_ind    = find(pairdata(:,1)==E2(1,i));
        if ~isempty(p_ind)
            xvn2(1,i) = xvn2(1,i) - pairdata(p_ind,2);
            yvn2(1,i) = yvn2(1,i) - pairdata(p_ind,3);
        end
        p_ind    = find(pairdata(:,1)==E3(1,i));
        if ~isempty(p_ind)
            xvn3(1,i) = xvn3(1,i) - pairdata(p_ind,2);
            yvn3(1,i) = yvn3(1,i) - pairdata(p_ind,3);
        end
    end
end

% compute face normals, lengths and mid-points
fnx = [yv2-yv1;yv3-yv2;yv1-yv3]; fny = -[xv2-xv1;xv3-xv2;xv1-xv3];
fL  = sqrt(fnx.^2 + fny.^2); 
fcx = [xv1+xv2; xv2+xv3; xv3+xv1]/2; fcy = [yv1+yv2; yv2+yv3; yv3+yv1]/2;

% Find boundary faces for each face 
id1 = find(BCType(:,1));
id2 = find(BCType(:,2));
id3 = find(BCType(:,3));

% Compute location of ghost neighbour vertices of reflected ghost elements at boundary faces
A0 = AVG2D*J*2;

H1 = 2*(A0(id1)./fL(1,id1)); % Using area of triangle formula
xvn1(id1) = xvn1(id1) + 2*(fnx(1,id1).*H1./fL(1,id1)); 
yvn1(id1) = yvn1(id1) + 2*(fny(1,id1).*H1./fL(1,id1)); 

H2 = 2*(A0(id2)./fL(2,id2)); 
xvn2(id2) = xvn2(id2) + 2*(fnx(2,id2).*H2./fL(2,id2));
yvn2(id2) = yvn2(id2) + 2*(fny(2,id2).*H2./fL(2,id2)); 

H3 = 2*(A0(id3)./fL(3,id3)); 
xvn3(id3) = xvn3(id3) + 2*(fnx(3,id3).*H3./fL(3,id3)); 
yvn3(id3) = yvn3(id3) + 2*(fny(3,id3).*H3./fL(3,id3)); 

% compute centroids for 4 triangles in each patch (note that the ghost
% nodes have been correctly evaluated at this point)
vcx  = [xv1+xv2+xv3; xv1+xv2+xvn1; xv2+xv3+xvn2; xv3+xv1+xvn3]/3;
vcy  = [yv1+yv2+yv3; yv1+yv2+yvn1; yv2+yv3+yvn2; yv3+yv1+yvn3]/3;


% Get the skew factor of original patch, and alphas needed by minmod
skew1        = zeros(3,K);
skew2        = zeros(3,K);
patch_alphas = zeros(3,3,K);
eps = 1e-12;
for i=1:K
    % Vector of centroid of central triangle to face mid-points  
    cfx = fcx(:,i) - vcx(1,i); cfy = fcy(:,i) - vcy(1,i); 
    cfL = sqrt(cfx.^2 + cfy.^2); 
    %cfx = cfx./cfL; cfy = cfy./cfL;
    
    % Vector of centroid of central triangle to centroids of nbs
    c2cx = (vcx(2:4,i) - vcx(1,i)); c2cy = (vcy(2:4,i) - vcy(1,i));
    c2cL = sqrt(c2cx.^2 + c2cy.^2); 
    %c2cx = c2cx./c2cL; c2cy = c2cy./c2cL;
    
    % Calculating minmod alphas
    Kind = [1,2,3,1,2];
    for j=1:3
        A = [c2cx(Kind(j)) , c2cx(Kind(j+1)); c2cy(Kind(j)) , c2cy(Kind(j+1))];
        b = [cfx(Kind(j)) ; cfy(Kind(j))];
        alphas = A\b;
        if(alphas(1) < -eps || alphas(2) < -eps)
            A = [c2cx(Kind(j)) , c2cx(Kind(j+2)); c2cy(Kind(j)) , c2cy(Kind(j+2))];
            alphas = A\b;
            assert(alphas(1) >= -eps & alphas(2) >= -eps)
            patch_alphas(j,:,i) = [alphas;Kind(j+2)];
        else
            patch_alphas(j,:,i) = [alphas;Kind(j+1)];
        end
    end

    % The two skew factors
    skew1(:,i) = (fnx(:,i)./fL(:,i)).*(cfx./cfL) + (fny(:,i)./fL(:,i)).*(cfy./cfL);
    skew2(:,i) = (fnx(:,i)./fL(:,i)).*(c2cx./c2cL) + (fny(:,i)./fL(:,i)).*(c2cy./c2cL); 
end



