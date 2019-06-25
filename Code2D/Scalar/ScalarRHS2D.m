function [rhsQ] = ScalarRHS2D(Q,time,mu,Problem,Mesh)

% function [rhsQ] = EulerRHS2D(Q,time, ExactSolutionBC);
% Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
%             with a local Lax-Friedrich flux

Mesh.vmapM = reshape(Mesh.vmapM, Mesh.Nfp*Mesh.Nfaces, Mesh.K); 
Mesh.vmapP = reshape(Mesh.vmapP, Mesh.Nfp*Mesh.Nfaces, Mesh.K);

% Solution traces
QP(:,:,1)=Q(Mesh.vmapP); QM(:,:,1)=Q(Mesh.vmapM);

% Apply boundary conditions to physical variable
QG = ApplyBCScalar2D(Q,time,Mesh,'phi');

% Set boundary conditions by modifying positive traces
if(~isempty(Mesh.mapBC_list))
        QPn = QP(:,:,1);
        QPn(Mesh.mapB) = QG(Mesh.Fmask(:,Mesh.Nfaces),:,1);
        QP(:,:,1) = QPn;
end

% Compute numerical fluxes of auxiliary equation
bd_flux_X=Mesh.LIFT*(Mesh.Fscale.*((QP+QM)/2.*Mesh.nx));
bd_flux_Y=Mesh.LIFT*(Mesh.Fscale.*((QP+QM)/2.*Mesh.ny));

%Compute internal fluxes for auxiliary variable
int_flux_X=Mesh.rx.*(Mesh.Drw*Q)+Mesh.sx.*(Mesh.Dsw*Q);
int_flux_Y=Mesh.ry.*(Mesh.Drw*Q)+Mesh.sy.*(Mesh.Dsw*Q);

% Compute auxiliary variables
qX=bd_flux_X-int_flux_X;
qY=bd_flux_Y-int_flux_Y;

% Add viscosity
gX=mu.*qX; gY=mu.*qY;

% Variable traces
gXP(:,:,1)=gX(Mesh.vmapP); gXM(:,:,1)=gX(Mesh.vmapM);
gYP(:,:,1)=gY(Mesh.vmapP); gYM(:,:,1)=gY(Mesh.vmapM);

% Apply boundary conditions to auxiliary variable
gXG = ApplyBCScalar2D(gX,time,Mesh,'aux');
gYG = ApplyBCScalar2D(gY,time,Mesh,'aux');

% Set boundary conditions by modifying positive traces
if(~isempty(Mesh.mapBC_list))
        gXPn = gXP(:,:,1);
        gXPn(Mesh.mapB) = gXG(Mesh.Fmask(:,Mesh.Nfaces),:,1);
        gXP(:,:,1) = gXPn;
        
        gYPn = gYP(:,:,1);
        gYPn(Mesh.mapB) = gYG(Mesh.Fmask(:,Mesh.Nfaces),:,1);
        gYP(:,:,1) = gYPn;
end

% Compute physical fluxes
[F,G] = ScalarFlux2D(Q, Problem);

% Compute weak derivatives
dFdr = Mesh.Drw*F(:,:,1); dFds = Mesh.Dsw*F(:,:,1);
dGdr = Mesh.Drw*G(:,:,1); dGds = Mesh.Dsw*G(:,:,1);

% Compute volume term of physical variable
int_flux = (Mesh.rx.*dFdr + Mesh.sx.*dFds) + (Mesh.ry.*dGdr + Mesh.sy.*dGds);
% Volume term for auxiliary variable
int_flux=int_flux-(Mesh.rx.*(Mesh.Drw*gX)+Mesh.sx.*(Mesh.Dsw*gX)+Mesh.ry.*(Mesh.Drw*gY)+Mesh.sy.*(Mesh.Dsw*gY));

% Evaluate primitive variables & flux functions 
[fM,gM] = ScalarFlux2D(QM, Problem);
[fP,gP] = ScalarFlux2D(QP, Problem);

% Compute maximum wave speed on the edges
lambda = max( Get_scalar_eig2D(QM,Problem),  ...
	          Get_scalar_eig2D(QP,Problem));
lambda = reshape(lambda, Mesh.Nfp, Mesh.Nfaces*Mesh.K);
lambda = ones(Mesh.Nfp, 1)*max(lambda, [], 1); 
lambda = reshape(lambda, Mesh.Nfp*Mesh.Nfaces, Mesh.K);

% Boundary flux for physical variable (Lax-Friedrichs)
flux = Mesh.nx.*(fP(:,:,1) + fM(:,:,1)) + Mesh.ny.*(gP(:,:,1) + gM(:,:,1)) + ...
       lambda.*(QM(:,:,1) - QP(:,:,1));
% Boundary flux from auxiliary variable
flux = flux - ((gXP+gXM).*Mesh.nx+(gYP+gYM).*Mesh.ny);
% Lift
bd_flux = Mesh.LIFT*(Mesh.Fscale.*flux/2);
   
% Compute rhs   
rhsQ(:,:,1) = int_flux - bd_flux;

return;

