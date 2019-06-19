function [rhsu] = ScalarRHS1D_weak(u,flux,dflux,mu,bc_cond,Mesh)

% Purpose  : Evaluate RHS flux in 1D Scalar equation using the weak form

% Apply BC to u
u_ext = Apply_BC1D(u(Mesh.VtoE),bc_cond);

% Compute numerical fluxes at interfaces for u
fluxr_u = ScalarC1D(u_ext(2,2:Mesh.K+1),u_ext(1,3:Mesh.K+2), @(u) u);
fluxl_u = ScalarC1D(u_ext(2,1:Mesh.K),u_ext(1,2:Mesh.K+1),@(u) u);

% Compute auxiliary variable variable q
qh = - Mesh.S'*u + (Mesh.Imat(:,Mesh.N+1)*fluxr_u(1,:) - Mesh.Imat(:,1)*fluxl_u(1,:));
q = mu.*Mesh.int_metric.*(Mesh.invM*qh); 

% Apply BC to auxiliary variable
bc_cond_aux=bc_cond;
bc_cond_aux{2}=0; bc_cond_aux{4}=0;
if strcmp(bc_cond{1},'N') 
    bc_cond_aux{1}='D';
else
    bc_cond_aux{1}='N';
end
if strcmp(bc_cond{3},'N') 
    bc_cond_aux{3}='D';
else
    bc_cond_aux{3}='N';
end 
q_ext = Apply_BC1D(q(Mesh.VtoE),bc_cond_aux);

% Compute numerical fluxes at interfaces for for mu*q and f
fluxr_q = ScalarC1D(q_ext(2,2:Mesh.K+1),q_ext(1,3:Mesh.K+2), @(u) u);
fluxl_q = ScalarC1D(q_ext(2,1:Mesh.K),q_ext(1,2:Mesh.K+1), @(u) u);
fluxr_fu = ScalarLF1D(u_ext(2,2:Mesh.K+1),u_ext(1,3:Mesh.K+2),flux,dflux);
fluxl_fu = ScalarLF1D(u_ext(2,1:Mesh.K),u_ext(1,2:Mesh.K+1),flux,dflux);

% compute right hand sides of the semi-discrete PDE
rhsu  = Mesh.S'*flux(u) - Mesh.S'*q + (Mesh.Imat(:,Mesh.Np)*fluxr_q(1,:) - Mesh.Imat(:,1)*fluxl_q(1,:)) - (Mesh.Imat(:,Mesh.Np)*fluxr_fu(1,:) - Mesh.Imat(:,1)*fluxl_fu(1,:)) ;
rhsu  = Mesh.int_metric.*(Mesh.invM*rhsu);

return