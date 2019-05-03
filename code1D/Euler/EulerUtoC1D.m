function [qc,Char] = EulerUtoC1D(q,gas_gamma,gas_const,Mesh)

% function to converts conserved variables to the characeristic variables
% for the Euler equations, using the cell average values of the 
% conserved variables

[NTC,dummy] = size(q(1,:,1)'); % All cells  on which q is defined. This may contain 
                      % ghost cells                     

qc = zeros(NTC,3); Char = zeros(Mesh.Np,NTC,3);

% Compute cell averages
% rhoh  = invV*q(:,:,1); rhoh(2:Np,:) =0; rhoa  = V*rhoh;  qc(:,1)  = rhoa(1,:);
% rhouh = invV*q(:,:,2); rhouh(2:Np,:)=0; rhoua = V*rhouh; qc(:,2)  = rhoua(1,:);
% Enerh = invV*q(:,:,3); Enerh(2:Np,:)=0; Enera = V*Enerh; qc(:,3)  = Enera(1,:);
qc(:,1) = Mesh.AVG1D*q(:,:,1);
qc(:,2) = Mesh.AVG1D*q(:,:,2);
qc(:,3) = Mesh.AVG1D*q(:,:,3);

% Compute characterisic variables
for i=1:NTC
    [L,invL] = EulerCharMat1D(qc(i,1),qc(i,2),qc(i,3),gas_gamma,gas_const);
    Char(:,i,:) = [q(:,i,1) q(:,i,2) q(:,i,3)]*invL';
end

return
