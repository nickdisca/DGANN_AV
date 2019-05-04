% Purpose : Compute outward pointing normals at elements faces and surface Jacobians

xr = Mesh.Dr*Mesh.x; yr = Mesh.Dr*Mesh.y; xs = Mesh.Ds*Mesh.x; ys = Mesh.Ds*Mesh.y; 
%Mesh.J  = xr.*ys-xs.*yr;

% interpolate geometric factors to face nodes
fxr = xr(Mesh.Fmask, :); fxs = xs(Mesh.Fmask, :); fyr = yr(Mesh.Fmask, :); fys = ys(Mesh.Fmask, :);

% build normals
Mesh.nx = zeros(3*Mesh.Nfp, Mesh.K); Mesh.ny = zeros(3*Mesh.Nfp, Mesh.K);
fid1 = (1:Mesh.Nfp)'; fid2 = (Mesh.Nfp+1:2*Mesh.Nfp)'; fid3 = (2*Mesh.Nfp+1:3*Mesh.Nfp)';

% face 1
Mesh.nx(fid1, :) =  fyr(fid1, :); Mesh.ny(fid1, :) = -fxr(fid1, :);

% face 2
Mesh.nx(fid2, :) =  fys(fid2, :)-fyr(fid2, :); Mesh.ny(fid2, :) = -fxs(fid2, :)+fxr(fid2, :);

% face 3
Mesh.nx(fid3, :) = -fys(fid3, :); Mesh.ny(fid3, :) =  fxs(fid3, :);

% normalise
Mesh.sJ = sqrt(Mesh.nx.^2+Mesh.ny.^2); Mesh.nx = Mesh.nx./Mesh.sJ; Mesh.ny = Mesh.ny./Mesh.sJ;

