function CFL = EulerCFL2D(Q, dt, Problem, Mesh)

% function dt = EulerCFL2D(Q, gas_gamma)
% purpose: compute the effective CFL corresponding to time step
% for the compressible Euler equations

rho = Q(:,:,1); rhou = Q(:,:,2); rhov = Q(:,:,3); Ener = Q(:,:,4);
rho = rho(Mesh.vmapM); rhou = rhou(Mesh.vmapM); 
rhov = rhov(Mesh.vmapM); Ener = Ener(Mesh.vmapM);

u = rhou./rho; v = rhov./rho;
p = (Problem.gas_gamma-1.0)*(Ener - rho.*(u.^2+v.^2)/2); c = sqrt(abs(Problem.gas_gamma*p./rho));

CFL = dt*max( ((Mesh.N+1)^2)*.5*Mesh.Fscale(:).*(sqrt ( u(:).^2 + v(:).^2 ) + c(:)));

return
