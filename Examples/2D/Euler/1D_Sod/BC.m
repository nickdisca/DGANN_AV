function [rhoG,uG,vG,preG] = BC(ckey,time,gas_gamma,gas_const)
 
% Purpose: Set Boundary Condition

p_L = 1.0;  rho_L = 1.0;    u_L = 0.0;     v_L = 0;
p_R = 0.1;  rho_R = 0.125;  u_R = 0.0;     v_R = 0;

if(ckey == 101)
    rhoG =@(x,y) rho_L*ones(size(x));
    uG   =@(x,y) u_L*ones(size(x));
    vG   =@(x,y) v_L*ones(size(x));
    preG =@(x,y) p_L*ones(size(x));
elseif(ckey == 102)
    rhoG =@(x,y) rho_R*ones(size(x));
    uG   =@(x,y) u_R*ones(size(x));
    vG   =@(x,y) v_R*ones(size(x));
    preG =@(x,y) p_R*ones(size(x));
end

return;


