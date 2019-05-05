% Find neighbors in patch
E1 = Mesh.EToE(:,1)'; E2 = Mesh.EToE(:,2)'; E3 = Mesh.EToE(:,3)';

% Find local faces in neighbouring triangles
f1 = Mesh.EToF(:,1)'; f2 = Mesh.EToF(:,2)'; f3 = Mesh.EToF(:,3)';

% Extract coordinates of vertices and centers of elements
v1 = Mesh.EToV(:,1); xv1 = Mesh.VX(v1); yv1 = Mesh.VY(v1);
v2 = Mesh.EToV(:,2); xv2 = Mesh.VX(v2); yv2 = Mesh.VY(v2);
v3 = Mesh.EToV(:,3); xv3 = Mesh.VX(v3); yv3 = Mesh.VY(v3);  


% Get vertices opposite shared face in neighbour
vind  = [3,1,2];
xvn1  = zeros(1,Mesh.K); xvn2 = xvn1; xvn3 = xvn1;
yvn1  = zeros(1,Mesh.K); yvn2 = yvn1; yvn3 = yvn1;
for i = 1:Mesh.K
   vn1   = Mesh.EToV(E1(1,i),vind(f1(1,i))); xvn1(1,i) = Mesh.VX(vn1); yvn1(1,i) = Mesh.VY(vn1);
   vn2   = Mesh.EToV(E2(1,i),vind(f2(1,i))); xvn2(1,i) = Mesh.VX(vn2); yvn2(1,i) = Mesh.VY(vn2);
   vn3   = Mesh.EToV(E3(1,i),vind(f3(1,i))); xvn3(1,i) = Mesh.VX(vn3); yvn3(1,i) = Mesh.VY(vn3);
   
   % Correctional shift periodic ghost triangles
   if(isKey(Mesh.PShift,i) && Mesh.UseMeshPerData)
        pairdata = Mesh.PShift(i);
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
id1 = find(Mesh.BCTag(:,1));
id2 = find(Mesh.BCTag(:,2));
id3 = find(Mesh.BCTag(:,3));

% Compute location of ghost neighbour vertices of reflected ghost elements at boundary faces
A0 = Mesh.AVG2D*Mesh.J*2;

H1 = 2*(A0(id1)./fL(1,id1)); % Using area of triangle formula
xvn1(id1) = xvn1(id1) + 2*(fnx(1,id1).*H1./fL(1,id1)); 
yvn1(id1) = yvn1(id1) + 2*(fny(1,id1).*H1./fL(1,id1)); 

H2 = 2*(A0(id2)./fL(2,id2)); 
xvn2(id2) = xvn2(id2) + 2*(fnx(2,id2).*H2./fL(2,id2));
yvn2(id2) = yvn2(id2) + 2*(fny(2,id2).*H2./fL(2,id2)); 

H3 = 2*(A0(id3)./fL(3,id3)); 
xvn3(id3) = xvn3(id3) + 2*(fnx(3,id3).*H3./fL(3,id3)); 
yvn3(id3) = yvn3(id3) + 2*(fny(3,id3).*H3./fL(3,id3)); 

% Generating ghost elements
Mesh.EToGE = zeros(size(Mesh.EToE));
Mesh.KG    = length(id1) + length(id2) + length(id3);
Mesh.xG = zeros(Mesh.Np,Mesh.KG);
Mesh.yG = zeros(Mesh.Np,Mesh.KG);

Mesh.EToGE(id1,1) = (1:length(id1))';
Mesh.xG(:,1:length(id1)) = 0.5*(-(Mesh.r+Mesh.s)*xv1(id1)+(1+Mesh.r)*xvn1(id1)+(1+Mesh.s)*xv2(id1));
Mesh.yG(:,1:length(id1)) = 0.5*(-(Mesh.r+Mesh.s)*yv1(id1)+(1+Mesh.r)*yvn1(id1)+(1+Mesh.s)*yv2(id1));

Mesh.EToGE(id2,2) = (length(id1)+1:length(id1)+length(id2))';
Mesh.xG(:,length(id1)+1:length(id1)+length(id2)) = 0.5*(-(Mesh.r+Mesh.s)*xv2(id2)+(1+Mesh.r)*xvn2(id2)+(1+Mesh.s)*xv3(id2));
Mesh.yG(:,length(id1)+1:length(id1)+length(id2)) = 0.5*(-(Mesh.r+Mesh.s)*yv2(id2)+(1+Mesh.r)*yvn2(id2)+(1+Mesh.s)*yv3(id2));

Mesh.EToGE(id3,3) = (length(id1)+length(id2)+1:length(id1)+length(id2)+length(id3))';
Mesh.xG(:,length(id1)+length(id2)+1:length(id1)+length(id2)+length(id3)) = ...
               0.5*(-(Mesh.r+Mesh.s)*xv3(id3)+(1+Mesh.r)*xvn3(id3)+(1+Mesh.s)*xv1(id3));
Mesh.yG(:,length(id1)+length(id2)+1:length(id1)+length(id2)+length(id3)) = ...
               0.5*(-(Mesh.r+Mesh.s)*yv3(id3)+(1+Mesh.r)*yvn3(id3)+(1+Mesh.s)*yv1(id3));


% compute centroids for 4 triangles in each patch (note that the ghost
% nodes have been correctly evaluated at this point)
vcx  = [xv1+xv2+xv3; xv1+xv2+xvn1; xv2+xv3+xvn2; xv3+xv1+xvn3]/3;
vcy  = [yv1+yv2+yv3; yv1+yv2+yvn1; yv2+yv3+yvn2; yv3+yv1+yvn3]/3;


% Get the skew factor of original patch, and alphas needed by minmod
% skew1        = zeros(3,K);
% skew2        = zeros(3,K);
Mesh.patch_alphas = zeros(3,3,Mesh.K);
eps = 1e-12;
for i=1:Mesh.K
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
%             if(~(alphas(1) >= -eps && alphas(2) >= -eps))
%                 A
%                 b
%                 alphas
%             end
            A = [c2cx(Kind(j)) , c2cx(Kind(j+2)); c2cy(Kind(j)) , c2cy(Kind(j+2))];
            alphas = A\b;
            assert((alphas(1) >= -eps & alphas(2) >= -eps),...
                'The mesh traingulation is not appropriate');
            Mesh.patch_alphas(j,:,i) = [alphas;Kind(j+2)];
        else
            Mesh.patch_alphas(j,:,i) = [alphas;Kind(j+1)];
        end
    end

    % The two skew factors
%     skew1(:,i) = (fnx(:,i)./fL(:,i)).*(cfx./cfL) + (fny(:,i)./fL(:,i)).*(cfy./cfL);
%     skew2(:,i) = (fnx(:,i)./fL(:,i)).*(c2cx./c2cL) + (fny(:,i)./fL(:,i)).*(c2cy./c2cL); 
end

% Creating mirroring indices mapping for ghost elements using a dummy
% triangle
Mesh.MMAP = zeros(Mesh.Np,3);
xx = [-1,1,-1]; yy=[-1,-1,1];
vindm = [1,3,2;2,1,3;3,2,1];
xo = 0.5*(-(Mesh.r+Mesh.s)*xx(1)+(1+Mesh.r)*xx(2)+(1+Mesh.s)*xx(3));
yo = 0.5*(-(Mesh.r+Mesh.s)*yy(1)+(1+Mesh.r)*yy(2)+(1+Mesh.s)*yy(3));

for i=1:3
    xm = 0.5*(-(Mesh.r+Mesh.s)*xx(vindm(i,1))+(1+Mesh.r)*xx(vindm(i,2))+(1+Mesh.s)*xx(vindm(i,3)));
    ym = 0.5*(-(Mesh.r+Mesh.s)*yy(vindm(i,1))+(1+Mesh.r)*yy(vindm(i,2))+(1+Mesh.s)*yy(vindm(i,3)));
    D = (xo -xm').^2 + (yo-ym').^2;
    [idM, idP] = find(sqrt(abs(D))<eps);
    Mesh.MMAP(:,i) = idM;
end




