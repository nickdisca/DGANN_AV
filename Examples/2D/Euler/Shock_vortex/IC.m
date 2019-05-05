function Q = IC(x, y, gas_gamma,gas_const)
 
% function Q = IC(x, y)
% Purpose: Set Initial condition for 2D Euler equations. This corresponds
% to left moving shock and right moving vortex


rho_L = 1.0;
u_L   = sqrt(gas_gamma);
v_L   = 0;
p_L   = 1.0;

p_R   = 1.3;
rho_R = rho_L*((gas_gamma+1)*p_R + gas_gamma-1)/((gas_gamma-1)*p_R + gas_gamma+1);
u_R   = sqrt(gas_gamma) + sqrt(2)*(1-p_R)/sqrt(gas_gamma - 1 + p_R*(gas_gamma+1));
v_R   = 0.0;

eps   = 0.3;
rc    = 0.05;
beta  = 0.204;
xc    = 0.25;
yc    = 0.5;
r2    = ((x-xc).^2 + (y-yc).^2)/rc^2;
phir  = exp(beta*(1-r2));
delu  = eps*(y-yc)/rc.*phir;
delv  = -eps*(x-xc)/rc.*phir;
delT  = -(gas_gamma-1)*eps^2/(4*beta*gas_gamma)*phir.^2;

% min(min(delu))
% max(max(delu))
% min(min(delv))
% max(max(delv))
% min(min(delT))
% max(max(delT))

T_L_mod   = p_L/gas_const/rho_L + delT;
rho_L_mod = (gas_const*T_L_mod).^(1/(gas_gamma-1));
p_L_mod   = rho_L_mod.^gas_gamma;

u     = (u_L + delu).*(x<0.5) + u_R*(x>=0.5); 
v     = (v_L + delv).*(x<0.5) + v_R*(x>=0.5);
rho   = (rho_L_mod).*(x<0.5) + rho_R*(x>=0.5);
pre   = (p_L_mod).*(x<0.5) + p_R*(x>=0.5);

Q(:,:,1) = rho;
Q(:,:,2) = rho.*u;
Q(:,:,3) = rho.*v;
Q(:,:,4) = Euler_Energy2D(rho,u,v,pre,gas_gamma);



return;


