function [rhsu] = ScalarRHS1D_weak(u,flux,dflux)

% Purpose  : Evaluate RHS flux in 1D Scalar equation using the weak form

Globals1D_DG;

u_ext = apply_bc(u(VtoE));

% Compute numerical fluxes at interfaces
fluxr = ScalarLF(u_ext(2,2:K+1),u_ext(1,3:K+2),flux,dflux);
fluxl = ScalarLF(u_ext(2,1:K),u_ext(1,2:K+1),flux,dflux);

% compute right hand sides of the semi-discrete PDE
rhsu  = S'*flux(u) - (Imat(:,Np)*fluxr(1,:) - Imat(:,1)*fluxl(1,:)) ;
rhsu  = int_metric.*(invM*rhsu);
return
