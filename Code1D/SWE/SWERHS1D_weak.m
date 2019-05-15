function rhsq = SWERHS1D_weak(q, gravity,bc_cond,Mesh)

% Purpose  : Evaluate RHS flux in 1D Shallow water with weak form

%Globals1D_DG;

depth = q(:,:,1); discharge = q(:,:,2);

% Extebd face data
depth_ext     =  Apply_BC1D(depth(Mesh.VtoE),bc_cond(1,:));
discharge_ext =  Apply_BC1D(discharge(Mesh.VtoE),bc_cond(2,:));

% Compute volume fluxes
fdep = discharge; fdis = 0.5*gravity*depth.^2 + discharge.^2./depth;

% Compute surface fluxes
fluxr = SWE_LF1D([depth_ext(2,2:Mesh.K+1)' discharge_ext(2,2:Mesh.K+1)'],...
               [depth_ext(1,3:Mesh.K+2)' discharge_ext(1,3:Mesh.K+2)'],...
               gravity)';
fdepr = fluxr(1,:); fdisr = fluxr(2,:);

fluxl = SWE_LF1D([depth_ext(2,1:Mesh.K)' discharge_ext(2,1:Mesh.K)'],...
               [depth_ext(1,2:Mesh.K+1)' discharge_ext(1,2:Mesh.K+1)'],...
               gravity)';           
fdepl = fluxl(1,:); fdisl = fluxl(2,:);

% Compute rhs
rhsu         = Mesh.S'*fdep - (Mesh.Imat(:,Mesh.Np)*fdepr(1,:) - Mesh.Imat(:,1)*fdepl(1,:)) ;
rhsq(:,:,1)  = Mesh.int_metric.*(Mesh.invM*rhsu);
rhsu         = Mesh.S'*fdis - (Mesh.Imat(:,Mesh.Np)*fdisr(1,:) - Mesh.Imat(:,1)*fdisl(1,:)) ;
rhsq(:,:,2)  = Mesh.int_metric.*(Mesh.invM*rhsu);

return
