function [qc,Char] = SWEUtoC1D(q,gravity,Mesh)

% function to converts conserved variables to the characeristic variables
% for the shallow water equations, using the cell average values of the 
% conserved variables

%Globals1D_DG;

NTC = size(q(1,:,1)); % All cells  on which q is defined. This may contain 
                      % ghost cells

qc = zeros(NTC,2); Char = zeros(Np,NTC,2); % Including ghost cells
% Compute cell averages
% depthh = invV*q(:,:,1); depthh(2:Np,:)=0; 
% deptha = V*depthh; qc(:,1) = deptha(1,:);
qc(:,1) = Mesh.AVG1D*q(:,:,1);

% dischargeh = invV*q(:,:,2); dischargeh(2:Np,:)=0; 
% dischargea = V*dischargeh; qc(:,2) = dischargea(1,:);
qc(:,2) = Mesh.AVG1D*q(:,:,2);


% Compute characterisic variables
for i=1:NTC
    [L,invL] = SWECharMat1D(qc(i,1),qc(i,2),gravity);
    Char(:,i,:) = [q(:,i,1) q(:,i,2)]*invL';
end

return;
