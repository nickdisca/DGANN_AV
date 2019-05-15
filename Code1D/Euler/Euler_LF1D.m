function fval = Euler_LF1D(u,v,gas_gamma, gas_const)

% Evauate local Lax Friedrich flux for 1D Euler equations
% u and v are conserved variables

% flux for u
pu   = (gas_gamma - 1)*(u(:,3) - 0.5*u(:,2).^2./u(:,1));
velu = u(:,2)./u(:,1);
fu   = [u(:,2) (pu + u(:,1).*velu.^2) (u(:,3) + pu).*velu ];

% flux for v
pv   = (gas_gamma - 1)*(v(:,3) - 0.5*v(:,2).^2./v(:,1));
velv = v(:,2)./v(:,1);
fv   = [v(:,2) (pv + v(:,1).*velv.^2) (v(:,3) + pv).*velv ];

% Maximum eigenvalue at face
lamu = sqrt(gas_gamma*pu./u(:,1)) + abs(velu);
lamv = sqrt(gas_gamma*pv./v(:,1)) + abs(velv);
lam  = max(lamu,lamv);

fval = 0.5*(fu + fv - (lam*ones(1,3)).*(v-u));


return
