function q = SWECtoU1D(Char,qc,gravity,Mesh)

% function to converts charactersitic variables to the conserved variables
% for the shallow water equations, using the cell average values of the 
% conserved variables

%Globals1D_DG;

NTC = size(qc(1,:,1)); % All cells  on which q is defined. This may contain 
                       % ghost cells

q = zeros(Mesh.Np,NTC,2);

% Compute conserved variables
for i=1:NTC
    [L,invL] = SWECharMat1D(qc(i,1),qc(i,2),gravity);
    Con = [Char(:,i,1) Char(:,i,2)]*L';
    q(:,i,1) = Con(:,1); q(:,i,2) = Con(:,2);
end

return;
