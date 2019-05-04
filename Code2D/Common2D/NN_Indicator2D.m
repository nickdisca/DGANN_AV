function [ind] = NN_Indicator2D(Q,QG,Mesh,Net,nconst_ind)


% Get Neighbours
E1 = Mesh.EToE(:,1)'; E2 = Mesh.EToE(:,2)'; E3 = Mesh.EToE(:,3)';

% Get modal coefficients of patch
Qm0 = Mesh.invV*Q; Qm1 = Mesh.invV*Q(:,E1); Qm2 = Mesh.invV*Q(:,E2); Qm3 = Mesh.invV*Q(:,E3);

% Replacing boundary element neighbours with ghost neighbours
GE1 = find(Mesh.EToGE(:,1))';
GE2 = find(Mesh.EToGE(:,2))';
GE3 = find(Mesh.EToGE(:,3))';
Qm1(:,GE1) = Mesh.invV*QG(:,Mesh.EToGE(GE1,1));
Qm2(:,GE2) = Mesh.invV*QG(:,Mesh.EToGE(GE2,2));
Qm3(:,GE3) = Mesh.invV*QG(:,Mesh.EToGE(GE3,3));

p1_ind = [1,2,Mesh.N+2];
X      = [Qm0(p1_ind,nconst_ind); Qm1(p1_ind,nconst_ind); Qm2(p1_ind,nconst_ind); Qm3(p1_ind,nconst_ind)];

bc_prob             = zeros(1,Mesh.K);
bc_prob(nconst_ind) = ind_MLP2D(X,Net);

ind = find(bc_prob > 0.5);


return;
