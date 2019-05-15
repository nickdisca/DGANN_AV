function fval = SWE_LF1D(u,v,gravity)

% Evauate local Lax Friedrich flux for 1D Shallow wtaer equations

% flux for u
fu = [u(:,2) 0.5*gravity*u(:,1).^2 + u(:,2).^2./u(:,1)];

% flux for v
fv = [v(:,2) 0.5*gravity*v(:,1).^2 + v(:,2).^2./v(:,1)];

% Maximum eigenvalue at face
lamu = sqrt(gravity*u(:,1)) + abs(u(:,2)./u(:,1));
lamv = sqrt(gravity*v(:,1)) + abs(v(:,2)./v(:,1));
lam  = max(lamu,lamv);

fval = 0.5*(fu + fv - (lam*ones(1,2)).*(v-u));


return
