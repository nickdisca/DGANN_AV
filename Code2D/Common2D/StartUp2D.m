% Purpose : Setup script, building operators, grid, metric, and connectivity tables.
% Definition of constants
Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;

% Compute nodal set
fprintf('... generating nodes\n')
[x,y] = Nodes2D(N); [r,s] = xytors(x,y);

% Build reference element matrices
fprintf('... generating basic matrices\n')
V          = Vandermonde2D(N,r,s); 
invV       = inv(V);
MassMatrix = invV'*invV;
[Dr,Ds]    = Dmatrices2D(N, r, s, V);

% build coordinates of all the nodes
fprintf('... generating nodes cordinates\n')
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
x = 0.5*(-(r+s)*VX(va)+(1+r)*VX(vb)+(1+s)*VX(vc));
y = 0.5*(-(r+s)*VY(va)+(1+r)*VY(vb)+(1+s)*VY(vc));

% find all the nodes that lie on each edge
fprintf('... generating cell face mask\n')
fmask1   = find( abs(s+1) < NODETOL)'; 
fmask2   = find( abs(r+s) < NODETOL)';
fmask3   = find( abs(r+1) < NODETOL)';
Fmask    = [fmask1;fmask2;fmask3]';
Fx = x(Fmask(:), :); Fy = y(Fmask(:), :);

% Create surface integral terms
fprintf('... generating face matrices\n')
[LIFT,M1D_1, M1D_2, M1D_3] = Lift2D();

% Creating averaging matrices
AVG2D   = sum(MassMatrix)/2;
AVG1D_1 = sum(M1D_1)/2; 
AVG1D_2 = sum(M1D_2)/2; 
AVG1D_3 = sum(M1D_3)/2;

% calculate geometric factors
fprintf('... calculating geometric transform factors\n')
[rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds);

% calculate geometric factors
fprintf('... generating face normal data\n')
[nx, ny, sJ] = Normals2D();
Fscale = sJ./(J(Fmask,:));

% calculate incircle radius for each triangle
fprintf('... calculating radius of incircles\n')
xscale2D;

% Build connectivity matrix
%[EToE, EToF] = tiConnect2D(EToV);
fprintf('... creating connectivity matrices\n')
[EToE, EToF, PShift] = Connect2D(EToV,BFaces,PerBToB_map,PerBFToF_map,...
                                 UseMeshPerData,VX,VY);

% Build connectivity maps
fprintf('... building face maps\n')
BuildMaps2D;

% Compute weak operators (could be done in preprocessing to save time)'
% NOTE that there is a transponse in the weak formulation, thus the
% following matrix form is obtained after multiplying by the inverse
% of the mass matrix
fprintf('... generating weak operators\n')
[Vr, Vs] = GradVandermonde2D(N, r, s);
Drw = (V*Vr')/(V*V'); Dsw = (V*Vs')/(V*V');

% Find projection matrices need for Fu-Shu indicator
%ProjectFromNb2D = Get_Projection_2D;
