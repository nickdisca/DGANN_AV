function [rhoB,uB,vB,preB] =  EulerOutflow2D(rhom,um,vm,prem,rhop,up,vp,prep, nx,ny,gas_gamma)

% Evaluate averaged eigenvalues
am  = sqrt(gas_gamma*prem./rhom);
ap  = sqrt(gas_gamma*prep./rhop);
a   = 0.5*(am+ap);
Unm = um.*nx + vm.*ny;
Unp = up.*nx + vp.*ny;
Un  = 0.5*(Unm + Unp);

% Supersonic outflow
rhoB = rhom; uB = um; vB = vm; preB = prem;

% Subsonic inflow (note that Un > 0)
ind = find(Un < a);
preB(ind) = prep(ind);
rhoB(ind) = rhom(ind) + (preB(ind) - prem(ind))./am(ind).^2;
uB(ind)   = um(ind) - nx(ind).*(preB(ind) - prep(ind))./am(ind)./rhom(ind);
vB(ind)   = vm(ind) - ny(ind).*(preB(ind) - prep(ind))./am(ind)./rhom(ind);

return