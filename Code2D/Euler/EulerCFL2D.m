function CFL = EulerCFL2D(Q, dt, gas_gamma)

% function dt = EulerCFL2D(Q, gas_gamma)
% purpose: compute the effective CFL corresponding to time step
% for the compressible Euler equations

Globals2D_DG;

rho = Q(:,:,1); rhou = Q(:,:,2); rhov = Q(:,:,3); Ener = Q(:,:,4);
rho = rho(vmapM); rhou = rhou(vmapM); rhov = rhov(vmapM); Ener = Ener(vmapM);

u = rhou./rho; v = rhov./rho;
p = (gas_gamma-1.0)*(Ener - rho.*(u.^2+v.^2)/2); c = sqrt(abs(gas_gamma*p./rho));

CFL = dt*max( ((N+1)^2)*.5*Fscale(:).*(sqrt ( u(:).^2 + v(:).^2 ) + c(:)));

%rhoprange = [min(min(rho)), max(max(rho)), min(min(p)), max(max(p))]
return
