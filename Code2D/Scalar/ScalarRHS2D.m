function [rhsQ] = ScalarRHS2D(Q,time,Problem,Mesh)

% function [rhsQ] = EulerRHS2D(Q,time, ExactSolutionBC);
% Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
%             with a local Lax-Friedrich flux

% Apply boundary conditions
QG = ApplyBCScalar2D(Q,time,Mesh);

Mesh.vmapM = reshape(Mesh.vmapM, Mesh.Nfp*Mesh.Nfaces, Mesh.K); 
Mesh.vmapP = reshape(Mesh.vmapP, Mesh.Nfp*Mesh.Nfaces, Mesh.K);

% 1. Compute volume contributions (NOW INDEPENDENT OF SURFACE TERMS)
[F,G] = ScalarFlux2D(Q, Problem);

% Compute weak derivatives
dFdr = Mesh.Drw*F(:,:,1); dFds = Mesh.Dsw*F(:,:,1);
dGdr = Mesh.Drw*G(:,:,1); dGds = Mesh.Dsw*G(:,:,1);
rhsQ(:,:,1) = (Mesh.rx.*dFdr + Mesh.sx.*dFds) + (Mesh.ry.*dGdr + Mesh.sy.*dGds);

    
% 2. Compute surface contributions 
% 2.1 evaluate '-' and '+' traces of conservative variables
QM(:,:,1) = Q(Mesh.vmapM); QP(:,:,1) = Q(Mesh.vmapP);


% % 2.2 set boundary conditions by modifying positive traces
if(~isempty(Mesh.mapBC_list))
        QPn = QP(:,:,1);
        QPn(mapB) = QG(Mesh.Fmask(:,Mesh.Nfaces),:,1);
        QP(:,:,1) = QPn;
end

% 2.3 evaluate primitive variables & flux functions at '-' and '+' traces
[fM,gM] = ScalarFlux2D(QM, Problem);
[fP,gP] = ScalarFlux2D(QP, Problem);

% 2.4 Compute local Lax-Friedrichs/Rusonov numerical fluxes
lambda = max( Get_scalar_eig2D(QM,Problem),  ...
	          Get_scalar_eig2D(QP,Problem));
lambda = reshape(lambda, Mesh.Nfp, Mesh.Nfaces*Mesh.K);
lambda = ones(Mesh.Nfp, 1)*max(lambda, [], 1); 
lambda = reshape(lambda, Mesh.Nfp*Mesh.Nfaces, Mesh.K);

% 2.5 Lift fluxes
flux = Mesh.nx.*(fP(:,:,1) + fM(:,:,1)) + Mesh.ny.*(gP(:,:,1) + gM(:,:,1)) + ...
       lambda.*(QM(:,:,1) - QP(:,:,1));
rhsQ(:,:,1) = rhsQ(:,:,1) - Mesh.LIFT*(Mesh.Fscale.*flux/2);

return;

