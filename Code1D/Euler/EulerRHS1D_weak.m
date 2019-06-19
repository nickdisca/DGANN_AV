function rhsq = EulerRHS1D_weak(aux, gas_gamma, gas_const,mu,bc_cond,Mesh)

% Purpose  : Evaluate RHS flux in 1D Euler with weak form

rho = aux(:,:,1); mmt = aux(:,:,2); Ener = aux(:,:,3); 
vel = mmt./rho; pre = (gas_gamma - 1)*(Ener - 0.5*rho.*vel.^2);

% Extend face data
rho_ext  =  Apply_BC1D(rho(Mesh.VtoE),bc_cond(1,:));
mmt_ext  =  Apply_BC1D(mmt(Mesh.VtoE),bc_cond(2,:));
Ener_ext =  Apply_BC1D(Ener(Mesh.VtoE),bc_cond(3,:));

% Compute numerical fluxes at interfaces for physical variables
fluxr_var = Euler_C1D([rho_ext(2,2:Mesh.K+1)' mmt_ext(2,2:Mesh.K+1)' Ener_ext(2,2:Mesh.K+1)'],...
                 [rho_ext(1,3:Mesh.K+2)' mmt_ext(1,3:Mesh.K+2)' Ener_ext(1,3:Mesh.K+2)'])';
fluxl_var = Euler_C1D([rho_ext(2,1:Mesh.K)'   mmt_ext(2,1:Mesh.K)'   Ener_ext(2,1:Mesh.K)'],...
                 [rho_ext(1,2:Mesh.K+1)' mmt_ext(1,2:Mesh.K+1)' Ener_ext(1,2:Mesh.K+1)'])';
       
% Compute auxiliary variable
for i=1:3
    auxh = - Mesh.S'*aux(:,:,i) + (Mesh.Imat(:,Mesh.N+1)*fluxr_var(i,:) - Mesh.Imat(:,1)*fluxl_var(i,:));
    aux(:,:,i) = mu(:,:,i).*Mesh.int_metric.*(Mesh.invM*auxh); 
end

% Apply BC to auxiliary variable
bc_cond_aux=bc_cond;
for i=1:3
    bc_cond_aux{i,2}=0; bc_cond_aux{i,4}=0;
    if strcmp(bc_cond{i,1},'N')
        bc_cond_aux{i,1}='D';
    else
        bc_cond_aux{i,1}='N';
    end
    if strcmp(bc_cond{i,3},'N')
        bc_cond_aux{i,3}='D';
    else
        bc_cond_aux{i,3}='N';
    end
end
tmp=aux(:,:,1); aux_ext_1 = Apply_BC1D(tmp(Mesh.VtoE),bc_cond_aux(1,:));
tmp=aux(:,:,2); aux_ext_2 = Apply_BC1D(tmp(Mesh.VtoE),bc_cond_aux(2,:));
tmp=aux(:,:,3); aux_ext_3 = Apply_BC1D(tmp(Mesh.VtoE),bc_cond_aux(3,:));

% Compute numerical fluxes at interfaces for auxiliary variable
fluxr_aux = Euler_C1D([aux_ext_1(2,2:Mesh.K+1)' aux_ext_2(2,2:Mesh.K+1)' aux_ext_3(2,2:Mesh.K+1)'],...
                 [aux_ext_1(1,3:Mesh.K+2)' aux_ext_2(1,3:Mesh.K+2)' aux_ext_3(1,3:Mesh.K+2)'])';
fluxl_aux = Euler_C1D([aux_ext_1(2,1:Mesh.K)'   aux_ext_2(2,1:Mesh.K)'   aux_ext_3(2,1:Mesh.K)'],...
                 [aux_ext_1(1,2:Mesh.K+1)' aux_ext_2(1,2:Mesh.K+1)' aux_ext_3(1,2:Mesh.K+1)'])';


% Compute volume fluxes
frho = mmt;
fmmt = pre + vel.*mmt;
fener = (Ener + pre).*vel;

% Compute real surface fluxes
fluxr = Euler_LF1D([rho_ext(2,2:Mesh.K+1)' mmt_ext(2,2:Mesh.K+1)' Ener_ext(2,2:Mesh.K+1)'],...
                 [rho_ext(1,3:Mesh.K+2)' mmt_ext(1,3:Mesh.K+2)' Ener_ext(1,3:Mesh.K+2)'],...
                 gas_gamma, gas_const)';
frhor = fluxr(1,:); fmmtr = fluxr(2,:); fenerr = fluxr(3,:);

fluxl = Euler_LF1D([rho_ext(2,1:Mesh.K)'   mmt_ext(2,1:Mesh.K)'   Ener_ext(2,1:Mesh.K)'],...
                 [rho_ext(1,2:Mesh.K+1)' mmt_ext(1,2:Mesh.K+1)' Ener_ext(1,2:Mesh.K+1)'],...
                 gas_gamma, gas_const)';

frhol = fluxl(1,:); fmmtl = fluxl(2,:); fenerl = fluxl(3,:);


% Compute rhs
rhsu         = Mesh.S'*frho - Mesh.S'*aux(:,:,1) - (Mesh.Imat(:,Mesh.Np)*frhor(1,:) - Mesh.Imat(:,1)*frhol(1,:)) + (Mesh.Imat(:,Mesh.Np)*fluxr_aux(1,:) - Mesh.Imat(:,1)*fluxl_aux(1,:));
rhsq(:,:,1)  = Mesh.int_metric.*(Mesh.invM*rhsu);
rhsu         = Mesh.S'*fmmt - Mesh.S'*aux(:,:,2) - (Mesh.Imat(:,Mesh.Np)*fmmtr(1,:) - Mesh.Imat(:,1)*fmmtl(1,:)) + (Mesh.Imat(:,Mesh.Np)*fluxr_aux(2,:) - Mesh.Imat(:,1)*fluxl_aux(2,:));
rhsq(:,:,2)  = Mesh.int_metric.*(Mesh.invM*rhsu);
rhsu         = Mesh.S'*fener - Mesh.S'*aux(:,:,3) - (Mesh.Imat(:,Mesh.Np)*fenerr(1,:) - Mesh.Imat(:,1)*fenerl(1,:)) + (Mesh.Imat(:,Mesh.Np)*fluxr_aux(3,:) - Mesh.Imat(:,1)*fluxl_aux(3,:));
rhsq(:,:,3)  = Mesh.int_metric.*(Mesh.invM*rhsu);

return
