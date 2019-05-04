function dt = ScalarDT2D(lam,CFL,Mesh)

% function dt = ScalarDT2D(maxeig)
% purpose: compute the time step dt for the 2D scalar conervation law

lam = lam(Mesh.vmapM); 
dt  = 1/max( ((Mesh.N+1)^2)*.5*Mesh.Fscale(:).*lam(:));
dt  = CFL*dt;

return
