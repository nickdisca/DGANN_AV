function q = EulerCtoU1D(Char,qc,gas_gamma,gas_const)

% function to converts charactersitic variables to the conserved variables
% for the Conserved equations, using the cell average values of the 
% conserved variables

Globals1D_DG;

q = zeros(Np,K,3);

% Compute conserved variables
for i=1:K
    [L,invL] = EulerCharMat(qc(i,1),qc(i,2),qc(i,3),gas_gamma,gas_const);
    Con = [Char(:,i,1) Char(:,i,2) Char(:,i,3)]*L';
    q(:,i,1) = Con(:,1); q(:,i,2) = Con(:,2); q(:,i,3) = Con(:,3);
end

return;
