% Purpose : Setup script, building operators, maps
Mesh.Np = Mesh.N+1;

% Compute basic Legendre Gauss Lobatto grid
fprintf('... generating nodes\n')
Mesh.r = JacobiGL(0,0,N);

% Build reference element matrices
fprintf('... generating basic matrices\n')
Mesh.V     = Vandermonde1D(Mesh.N, Mesh.r); 
Mesh.invV  = inv(Mesh.V);
Mesh.Dr    = Dmatrix1D(Mesh.N, Mesh.r, Mesh.V);
Mesh.M     = inv(Mesh.V')/Mesh.V; % mass matrix corresponding to reference element
Mesh.invM  = inv(Mesh.M);
Mesh.S     = Mesh.M*Mesh.Dr;
Mesh.int_metric = 2*ones(Mesh.Np,1)*(1./Mesh.hK);
Mesh.Imat  = eye(Mesh.Np);
Mesh.AVG1D = sum(Mesh.M)/2; % Averaging array

% build coordinates of all the nodes
fprintf('... generating nodes cordinates\n')
Mesh.x = ones(Mesh.N+1,1)*Mesh.VX(1:Mesh.Nv-1) ...
        + 0.5*(Mesh.r+1)*(Mesh.VX(2:Mesh.Nv)-Mesh.VX(1:Mesh.Nv-1));

% Create face node extractor
fprintf('... creating connectivity matrices\n')
Mesh.VtoE = zeros(2,Mesh.K);
for j=1:Mesh.K
    Mesh.VtoE(1,j) = (j-1)*(Mesh.Np)+1;
    Mesh.VtoE(2,j) = j*(Mesh.Np);
end

% Build Projection maps (needed for Shu-Fu indicator)
%Get_Projection_1D;



