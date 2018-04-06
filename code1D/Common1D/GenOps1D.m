% Purpose : Setup script, building operators, maps
Np = N+1;

% Compute basic Legendre Gauss Lobatto grid
r = JacobiGL(0,0,N);

% Build reference element matrices
V     = Vandermonde1D(N, r); 
invV  = inv(V);
Dr    = Dmatrix1D(N, r, V);
M     = inv(V')/V; % mass matrix corresponding to reference element
invM  = inv(M);
S     = M*Dr;
int_metric = 2*ones(Np,1)*(1./hK);
Imat  = eye(Np);

% Get extrapolations needed for Fu Shu indicator
% [extrap_left1D, extrap_right1D] = ExtrapolateBasis1D(N,r,2);

% build coordinates of all the nodes
x = ones(N+1,1)*VX(1:Nv-1) + 0.5*(r+1)*(VX(2:Nv)-VX(1:Nv-1));

% Create face node extractor
VtoE = zeros(2,K);
for j=1:K
    VtoE(1,j) = (j-1)*(Np)+1;
    VtoE(2,j) = j*(Np);
end



