% Purpose  : Compute surface to volume lift term for DG formulation

Emat = zeros(Mesh.Np, Mesh.Nfaces*Mesh.Nfp);
Vmid1D = Vandermonde1D(Mesh.N,0*Mesh.r);

% face 1
faceR = Mesh.r(Mesh.Fmask(:,1));
V1D = Vandermonde1D(Mesh.N, faceR); 
Mesh.M1D_1 = inv(V1D*V1D');
Emat(Mesh.Fmask(:,1),1:Mesh.Nfp) = Mesh.M1D_1;
Mesh.facemid1 = Vmid1D(1,:)/V1D;

% face 2
faceR = Mesh.r(Mesh.Fmask(:,2));
V1D = Vandermonde1D(Mesh.N, faceR);
Mesh.M1D_2 = inv(V1D*V1D');
Emat(Mesh.Fmask(:,2),Mesh.Nfp+1:2*Mesh.Nfp) = Mesh.M1D_2;
Mesh.facemid2 = Vmid1D(1,:)/V1D;

% face 3
faceS = Mesh.s(Mesh.Fmask(:,3));
V1D = Vandermonde1D(Mesh.N, faceS); 
Mesh.M1D_3 = inv(V1D*V1D');
Emat(Mesh.Fmask(:,3),2*Mesh.Nfp+1:3*Mesh.Nfp) = Mesh.M1D_3;
Mesh.facemid3 = Vmid1D(1,:)/V1D;

% inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
Mesh.LIFT = Mesh.V*(Mesh.V'*Emat);

%LIFT,M1D_1, M1D_2, M1D_3, facemid1, facemid2, facemid3