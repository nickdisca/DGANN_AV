function dt = ScalarDT2D(lam,visc,CFL,Mesh)

% function dt = ScalarDT2D(maxeig,maxvisc)
% purpose: compute the time step dt for the 2D scalar conervation law

dt = 1/(max(lam(:))*(Mesh.N^2)/min(2*Mesh.dx(:))+visc*(Mesh.N^4)/(min(2*Mesh.dx)^2));
dt = CFL*dt;

return
