function Q = IC(x, y, gas_gamma, gas_const)
 
% function Q = IC(x, y)
% Purpose: Set Initial conditions for double mach reflection
% Location of initial discontinuity
xc = 1/6; yc = 0; aos = pi/3;

p1 = 116.5;  rho1 = 8.0;    u1 = 7.14471;     v1 = -4.125;
p2 = 1.0;    rho2 = 1.4;    u2 = 0.0;     v2 = 0;

% Initial profile
pre = p1*(y>(x-xc)*tan(aos) + yc) + p2*(y<=(x-xc)*tan(aos) + yc);
u   = u1*(y>(x-xc)*tan(aos) + yc) + u2*(y<=(x-xc)*tan(aos) + yc);
v   = v1*(y>(x-xc)*tan(aos) + yc) + v2*(y<=(x-xc)*tan(aos) + yc);
rho = rho1*(y>(x-xc)*tan(aos) + yc) + rho2*(y<=(x-xc)*tan(aos) + yc);

Q(:,:,1) = rho; Q(:,:,2) = rho.*u; Q(:,:,3) = rho.*v;
Q(:,:,4) = Euler_Energy2D(rho,u,v,pre,gas_gamma);

return;


