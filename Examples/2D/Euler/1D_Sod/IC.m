function Q = IC(x, y, gas_gamma, gas_const)
 
% function Q = IC(x, y)
% Purpose: Set Initial conditio 
% Location of initial discontinuity
xc = 0; 

p1 = 1.0;  rho1 = 1.0;    u1 = 0.0;     v1 = 0;
p2 = 0.1;  rho2 = 0.125;  u2 = 0.0;     v2 = 0;

% Initial profile
pre = p1*(x<xc) + p2*(x>=xc);
u   = u1*(x<xc) + u2*(x>=xc);
v   = v1*(x<xc) + v2*(x>=xc);
rho = rho1*(x<xc) + rho2*(x>=xc);

Q(:,:,1) = rho; Q(:,:,2) = rho.*u; Q(:,:,3) = rho.*v;
Q(:,:,4) = Euler_Energy2D(rho,u,v,pre,gas_gamma);

return;


