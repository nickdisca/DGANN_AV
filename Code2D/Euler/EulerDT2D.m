function dt = EulerDT2D(Q, gas_gamma)

% function dt = EulerDT2D(Q, gas_gamma)
% purpose: compute the time step dt for the compressible Euler equations

Globals2D_DG;

rho = Q(:,:,1); rhou = Q(:,:,2); rhov = Q(:,:,3); Ener = Q(:,:,4);
rho = rho(vmapM); rhou = rhou(vmapM); rhov = rhov(vmapM); Ener = Ener(vmapM);

u = rhou./rho; v = rhov./rho;
p = (gas_gamma-1.0)*(Ener - rho.*(u.^2+v.^2)/2); c = sqrt(abs(gas_gamma*p./rho));

dt = CFL/max( ((N+1)^2)*.5*Fscale(:).*(sqrt ( u(:).^2 + v(:).^2 ) + c(:)));

%rhoprange = [min(min(rho)), max(max(rho)), min(min(p)), max(max(p))]
return
