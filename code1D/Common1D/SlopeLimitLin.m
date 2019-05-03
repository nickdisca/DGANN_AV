function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1,Mesh)

% function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
% Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
%          (vm1,v0,vp1) are cell averages left, center, and right

%Globals1D_DG;

% Compute various geometric measures
ulimit = ul; h = xl(Mesh.Np,:)-xl(1,:); 
x0 = ones(Mesh.Np,1)*(xl(1,:) + h/2);

hN = ones(Mesh.Np,1)*h;

% Limit function
ux = (2./hN).*(Mesh.Dr*ul);

ulimit = ones(Mesh.Np,1)*v0+(xl-x0).*(ones(Mesh.Np,1)*minmod([ux(1,:);(vp1-v0)./h;(v0-vm1)./h]));
return
