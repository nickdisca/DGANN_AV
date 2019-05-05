function [rhoG,uG,vG,preG] = BC(ckey,time,gas_gamma,gas_const)

% function [rhoG,uG,vG,preG] = BC(ckey,time,gas_gamma)
% Purpose: Set Boundary Condition
rho_L = 1.0;
u_L   = sqrt(gas_gamma);
v_L   = 0;
p_L   = 1.0;

p_R   = 1.3;
rho_R = rho_L*((gas_gamma+1)*p_R + gas_gamma-1)/((gas_gamma-1)*p_R + gas_gamma+1);
u_R   = sqrt(gas_gamma) + sqrt(2)*(1-p_R)/sqrt(gas_gamma - 1 + p_R*(gas_gamma+1));
v_R   = 0.0;

if(ckey == 101)
    rhoG =@(x,y) rho_L*ones(size(x));
    uG   =@(x,y) u_L*ones(size(x));
    vG   =@(x,y) v_L*ones(size(x));
    preG =@(x,y) p_L*ones(size(x));
elseif(ckey == 102)
    rhoG =@(x,y) 0*x;
    uG   =@(x,y) 0*x;
    vG   =@(x,y) 0*x;
    preG =@(x,y) 0*x;
elseif(ckey == 104)
    rhoG =@(x,y) 0*x;
    uG   =@(x,y) 0*x;
    vG   =@(x,y) 0*x;
    preG =@(x,y) 0*x;
elseif(ckey == 103)
    rhoG =@(x,y) rho_R*ones(size(x));
    uG   =@(x,y) u_R*ones(size(x));
    vG   =@(x,y) v_R*ones(size(x));
    preG =@(x,y) p_R*ones(size(x));
end



return;


