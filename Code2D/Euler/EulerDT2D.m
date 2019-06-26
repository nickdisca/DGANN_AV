function dt = EulerDT2D(lam,visc,CFL,Mesh)

% function dt = EulerDT2D(Q, gas_gamma)
% purpose: compute the time step dt for the compressible Euler equations

dt = 1/(max(lam(:))*(Mesh.N^2)/min(2*Mesh.dx(:))+visc*(Mesh.N^4)/(min(2*Mesh.dx)^2));
dt = CFL*dt;


return
