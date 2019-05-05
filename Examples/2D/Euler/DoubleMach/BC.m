function [rhoG,uG,vG,preG] = BC(ckey,time,gas_gamma,gas_const)

% function [rhoG,uG,vG,preG] = BC(ckey,time,gas_gamma)
% Purpose: Set Boundary Condition
xc = 1/6; yc = 0; aos = pi/3; M = 10;
pre1 = 116.5;  rho1 = 8.0;    u1 = 7.14471; v1 = -4.125;
pre2 = 1.0;    rho2 = 1.4;    u2 = 0.0;     v2 = 0;

xs = time*M/sin(aos);

if(ckey == 101 || ckey == 103)
    rhoG =@(x,y) rho1*ones(size(x));
    uG   =@(x,y) u1*ones(size(x));
    vG   =@(x,y) v1*ones(size(x));
    preG =@(x,y) pre1*ones(size(x));
elseif(ckey == 102)
    rhoG =@(x,y) rho2*ones(size(x));
    uG   =@(x,y) u2*ones(size(x));
    vG   =@(x,y) v2*ones(size(x));
    preG =@(x,y) pre2*ones(size(x));
elseif(ckey == 104)
    rhoG =@(x,y) 0*x;
    uG   =@(x,y) 0*x;
    vG   =@(x,y) 0*x;
    preG =@(x,y) 0*x;
elseif(ckey == 105)
    rhoG =@(x,y) rho1*(y>(x-xc-xs)*tan(aos) + yc) + rho2*(y<=(x-xc-xs)*tan(aos) + yc);
    uG   =@(x,y) u1*(y>(x-xc-xs)*tan(aos) + yc) + u2*(y<=(x-xc-xs)*tan(aos) + yc);
    vG   =@(x,y) v1*(y>(x-xc-xs)*tan(aos) + yc) + v2*(y<=(x-xc-xs)*tan(aos) + yc);
    preG =@(x,y) pre1*(y>(x-xc-xs)*tan(aos) + yc) + pre2*(y<=(x-xc-xs)*tan(aos) + yc);
end



return;


