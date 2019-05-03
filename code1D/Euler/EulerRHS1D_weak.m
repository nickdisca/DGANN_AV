function rhsq = EulerRHS1D_weak(q, gas_gamma, gas_const,bc_cond,Mesh)

% Purpose  : Evaluate RHS flux in 1D Shallow water with weak form

rho = q(:,:,1); mmt = q(:,:,2); Ener = q(:,:,3); 
vel = mmt./rho; pre = (gas_gamma - 1)*(Ener - 0.5*rho.*vel.^2);

% Extebd face data
rho_ext  =  Apply_BC1D(rho(Mesh.VtoE),bc_cond(1,:));
mmt_ext  =  Apply_BC1D(mmt(Mesh.VtoE),bc_cond(2,:));
Ener_ext =  Apply_BC1D(Ener(Mesh.VtoE),bc_cond(3,:));

% Compute volume fluxes
frho = mmt;
fmmt = pre + vel.*mmt;
fener = (Ener + pre).*vel;

% Compute surface fluxes
fluxr = Euler_LF1D([rho_ext(2,2:Mesh.K+1)' mmt_ext(2,2:Mesh.K+1)' Ener_ext(2,2:Mesh.K+1)'],...
                 [rho_ext(1,3:Mesh.K+2)' mmt_ext(1,3:Mesh.K+2)' Ener_ext(1,3:Mesh.K+2)'],...
                 gas_gamma, gas_const)';
frhor = fluxr(1,:); fmmtr = fluxr(2,:); fenerr = fluxr(3,:);

fluxl = Euler_LF1D([rho_ext(2,1:Mesh.K)'   mmt_ext(2,1:Mesh.K)'   Ener_ext(2,1:Mesh.K)'],...
                 [rho_ext(1,2:Mesh.K+1)' mmt_ext(1,2:Mesh.K+1)' Ener_ext(1,2:Mesh.K+1)'],...
                 gas_gamma, gas_const)';

frhol = fluxl(1,:); fmmtl = fluxl(2,:); fenerl = fluxl(3,:);


% Compute rhs
rhsu         = Mesh.S'*frho - (Mesh.Imat(:,Mesh.Np)*frhor(1,:) - Mesh.Imat(:,1)*frhol(1,:)) ;
rhsq(:,:,1)  = Mesh.int_metric.*(Mesh.invM*rhsu);
rhsu         = Mesh.S'*fmmt - (Mesh.Imat(:,Mesh.Np)*fmmtr(1,:) - Mesh.Imat(:,1)*fmmtl(1,:)) ;
rhsq(:,:,2)  = Mesh.int_metric.*(Mesh.invM*rhsu);
rhsu         = Mesh.S'*fener - (Mesh.Imat(:,Mesh.Np)*fenerr(1,:) - Mesh.Imat(:,1)*fenerl(1,:)) ;
rhsq(:,:,3)  = Mesh.int_metric.*(Mesh.invM*rhsu);

return
