function Q = IC(x, y, gas_gamma, gas_const)
 
% function Q = IC(x, y)
% Purpose: Set Initial conditio for 2D Riemann problem for Euler equations
% as per "solution of Two-Dimensiona Riemann Problems for Gas Dynamics
% without Riemann Problem Solvers" by Kurganov and Tadmor (2000)

% Location of initial discontinuity
xc = 0; yc = 0; gamma = 1.4;

% Config 4
p1 = 0.4;  rho1 = 0.5313;    u1 = 0;      v1 = 0;
p2 = 1.0;  rho2 = 1.0;       u2 = 0.7276; v2 = 0;
p3 = 1.0;  rho3 = 0.8;       u3 = 0.0;    v3 = 0.0;
p4 = 1.0;  rho4 = 1.0;       u4 = 0;      v4 = 0.7276;

% Initial profile
pre = p1*(x>xc).*(y>yc)   + p2*(x<=xc).*(y>yc)   + p3*(x<=xc).*(y<=yc)    + p4*(x>xc).*(y<=yc);
u   = u1*(x>xc).*(y>yc)   + u2*(x<=xc).*(y>yc)   + u3*(x<=xc).*(y<=yc)    + u4*(x>xc).*(y<=yc);
v   = v1*(x>xc).*(y>yc)   + v2*(x<=xc).*(y>yc)   + v3*(x<=xc).*(y<=yc)    + v4*(x>xc).*(y<=yc);
rho = rho1*(x>xc).*(y>yc) + rho2*(x<=xc).*(y>yc) + rho3*(x<=xc).*(y<=yc)  + rho4*(x>xc).*(y<=yc);

Q(:,:,1) = rho; Q(:,:,2) = rho.*u; Q(:,:,3) = rho.*v;
Q(:,:,4) = Euler_Energy2D(rho,u,v,pre,gas_gamma);

return;


