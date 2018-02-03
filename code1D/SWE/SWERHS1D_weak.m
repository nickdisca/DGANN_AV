function rhsq = SWERHS1D_weak(q, gravity)

% Purpose  : Evaluate RHS flux in 1D Shallow water with weak form

Globals1D_DG;

depth = q(:,:,1); discharge = q(:,:,2);

% Extebd face data
depth_ext     =  apply_bc(depth(VtoE));
discharge_ext =  apply_bc(discharge(VtoE));

% Compute volume fluxes
fdep = discharge; fdis = 0.5*gravity*depth.^2 + discharge.^2./depth;

% Compute surface fluxes
fluxr = SWE_LF([depth_ext(2,2:K+1)' discharge_ext(2,2:K+1)'],...
               [depth_ext(1,3:K+2)' discharge_ext(1,3:K+2)'],...
               gravity)';
fdepr = fluxr(1,:); fdisr = fluxr(2,:);

fluxl = SWE_LF([depth_ext(2,1:K)' discharge_ext(2,1:K)'],...
               [depth_ext(1,2:K+1)' discharge_ext(1,2:K+1)'],...
               gravity)';           
fdepl = fluxl(1,:); fdisl = fluxl(2,:);

% Compute rhs
rhsu         = S'*fdep - (Imat(:,Np)*fdepr(1,:) - Imat(:,1)*fdepl(1,:)) ;
rhsq(:,:,1)  = int_metric.*(invM*rhsu);
rhsu         = S'*fdis - (Imat(:,Np)*fdisr(1,:) - Imat(:,1)*fdisl(1,:)) ;
rhsq(:,:,2)  = int_metric.*(invM*rhsu);

return
