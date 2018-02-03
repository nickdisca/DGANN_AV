function [qc,Char] = EulerUtoC1D(q,gas_gamma, gas_const)

% function to converts conserved variables to the characeristic variables
% for the Euler equations, using the cell average values of the 
% conserved variables

Globals1D_DG;

qc = zeros(K,3); Char = zeros(Np,K,3);

% Compute cell averages
rhoh  = invV*q(:,:,1); rhoh(2:Np,:) =0; rhoa  = V*rhoh;  qc(:,1)  = rhoa(1,:);
rhouh = invV*q(:,:,2); rhouh(2:Np,:)=0; rhoua = V*rhouh; qc(:,2)  = rhoua(1,:);
Enerh = invV*q(:,:,3); Enerh(2:Np,:)=0; Enera = V*Enerh; qc(:,3)  = Enera(1,:);

% Compute characterisic variables
for i=1:K
    [L,invL] = EulerCharMat(qc(i,1),qc(i,2),qc(i,3),gas_gamma,gas_const);
    Char(:,i,:) = [q(:,i,1) q(:,i,2) q(:,i,3)]*invS';
end

return
