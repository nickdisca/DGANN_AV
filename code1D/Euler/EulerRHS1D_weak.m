function rhsq = EulerRHS1D_weak(q, gas_gamma, gas_const)

% Purpose  : Evaluate RHS flux in 1D Shallow water with weak form

Globals1D_DG;

rho = q(:,:,1); mmt = q(:,:,2); Ener = q(:,:,3); 
vel = mmt./rho; pre = (gas_gamma - 1)*(Ener - 0.5*rho.*vel.^2);

% Extebd face data
rho_ext  =  Apply_BC1D(rho(VtoE),bc_cond(1,:));
mmt_ext  =  Apply_BC1D(mmt(VtoE),bc_cond(2,:));
Ener_ext =  Apply_BC1D(Ener(VtoE),bc_cond(3,:));

% Compute volume fluxes
frho = mmt;
fmmt = pre + vel.*mmt;
fener = (Ener + pre).*vel;

% Compute surface fluxes
fluxr = Euler_LF1D([rho_ext(2,2:K+1)' mmt_ext(2,2:K+1)' Ener_ext(2,2:K+1)'],...
                 [rho_ext(1,3:K+2)' mmt_ext(1,3:K+2)' Ener_ext(1,3:K+2)'],...
                 gas_gamma, gas_const)';
frhor = fluxr(1,:); fmmtr = fluxr(2,:); fenerr = fluxr(3,:);

fluxl = Euler_LF1D([rho_ext(2,1:K)'   mmt_ext(2,1:K)'   Ener_ext(2,1:K)'],...
                 [rho_ext(1,2:K+1)' mmt_ext(1,2:K+1)' Ener_ext(1,2:K+1)'],...
                 gas_gamma, gas_const)';

frhol = fluxl(1,:); fmmtl = fluxl(2,:); fenerl = fluxl(3,:);


% Compute rhs
rhsu         = S'*frho - (Imat(:,Np)*frhor(1,:) - Imat(:,1)*frhol(1,:)) ;
rhsq(:,:,1)  = int_metric.*(invM*rhsu);
rhsu         = S'*fmmt - (Imat(:,Np)*fmmtr(1,:) - Imat(:,1)*fmmtl(1,:)) ;
rhsq(:,:,2)  = int_metric.*(invM*rhsu);
rhsu         = S'*fener - (Imat(:,Np)*fenerr(1,:) - Imat(:,1)*fenerl(1,:)) ;
rhsq(:,:,3)  = int_metric.*(invM*rhsu);

return
