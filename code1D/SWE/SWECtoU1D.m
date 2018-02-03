function q = SWECtoU1D(Char,qc,gravity)

% function to converts charactersitic variables to the conserved variables
% for the shallow water equations, using the cell average values of the 
% conserved variables

Globals1D_DG;

q = zeros(Np,K,2);

% Compute conserved variables
for i=1:K
    [L,invL] = SWECharMat(qc(i,1),qc(i,2),gravity);
    Con = [Char(:,i,1) Char(:,i,2)]*L';
    q(:,i,1) = Con(:,1); q(:,i,2) = Con(:,2);
end

return;
