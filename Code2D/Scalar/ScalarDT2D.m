function dt = ScalarDT2D(lam)

% function dt = ScalarDT2D(maxeig)
% purpose: compute the time step dt for the 2D scalar conervation law

Globals2D_DG;

lam = lam(vmapM); 
dt  = 1/max( ((N+1)^2)*.5*Fscale(:).*lam(:));
dt  = CFL*dt;

return
