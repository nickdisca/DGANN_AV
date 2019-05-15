% Purpose  : Generate a equidistant grid with K elements. If mesh_pert > 0
% then the mesh is perturbed randomly

function [Nv, VX, hK] = MeshGen1D(xmin,xmax,K,mesh_pert)

Nv         = K+1; 
h          = (xmax - xmin)/K;
pert_scale = mesh_pert*h; 

% Generate node coordinates
VX = (1:Nv);
for i = 1:Nv
  VX(i) = (xmax-xmin)*(i-1)/(Nv-1) + xmin;
end

% Perturb interior nodes
rng(1);
pert       = (rand(1,Nv-2) - 0.5)*pert_scale;
VX(2:Nv-1) = VX(2:Nv-1) + pert;
hK         = VX(2:Nv) - VX(1:Nv-1);

return
