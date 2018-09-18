function [rhoB,uB,vB,preB] =  EulerInflow2D(rhom,um,vm,prem,rhop,up,vp,prep, nx,ny,gas_gamma)


% Evaluate averaged eigenvalues
am  = sqrt(gas_gamma*prem./rhom);
ap  = sqrt(gas_gamma*prep./rhop);
a   = 0.5*(am+ap);
Unm = um.*nx + vm.*ny;
Unp = up.*nx + vp.*ny;
Un  = 0.5*(Unm + Unp);

% Supersonic inflow
rhoB = rhop; uB = up; vB = vp; preB = prep;

% Subsonic inflow (note that Un < 0)
ind = find(-Un < a);
preB(ind) =0.5*(prep(ind) + prem(ind) - rhom(ind).*am(ind).*(...
                nx(ind).*(up(ind)-um(ind)) + ny(ind).*(vp(ind)-vm(ind)) ) );
rhoB(ind) = rhop(ind) + (preB(ind) - prep(ind))./am(ind).^2;
uB(ind)   = up(ind) + nx(ind).*(preB(ind) - prep(ind))./am(ind)./rhom(ind);
vB(ind)   = vp(ind) + ny(ind).*(preB(ind) - prep(ind))./am(ind)./rhom(ind);


return