function [rhsu] = ScalarRHS1D_weak(u,flux,dflux,bc_cond,Mesh)

% Purpose  : Evaluate RHS flux in 1D Scalar equation using the weak form

%Globals1D_DG;

u_ext = Apply_BC1D(u(Mesh.VtoE),bc_cond);

% Compute numerical fluxes at interfaces
fluxr = ScalarLF1D(u_ext(2,2:Mesh.K+1),u_ext(1,3:Mesh.K+2),flux,dflux);
fluxl = ScalarLF1D(u_ext(2,1:Mesh.K),u_ext(1,2:Mesh.K+1),flux,dflux);

% compute right hand sides of the semi-discrete PDE
rhsu  = Mesh.S'*flux(u) - (Mesh.Imat(:,Mesh.Np)*fluxr(1,:) - Mesh.Imat(:,1)*fluxl(1,:)) ;
rhsu  = Mesh.int_metric.*(Mesh.invM*rhsu);
return
