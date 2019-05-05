function Q = IC(x, y, gas_gamma,gas_const)
 
% function Q = IC(x, y)
% Purpose: Set Initial condition for 2D Euler equations. This corresponds
% to the isentropic vortex

beta  = 5;
alpha = pi/4;
M     = 0.1;
xc    = 0;
yc    = 0;
r2    = (x-xc).^2 + (y-yc).^2;
rho = (1 - beta^2*(gas_gamma-1)*exp(1-r2)/(8*pi*pi*gas_gamma)).^(1/(gas_gamma-1));
u   = M*cos(alpha) - beta*(y-yc).*exp(0.5*(1-r2))/(2*pi);
v   = M*sin(alpha) + beta*(x-xc).*exp(0.5*(1-r2))/(2*pi);
pre = rho.^gas_gamma;

Q(:,:,1) = rho;
Q(:,:,2) = rho.*u;
Q(:,:,3) = rho.*v;
Q(:,:,4) = Euler_Energy2D(rho,u,v,pre,gas_gamma);
return;


