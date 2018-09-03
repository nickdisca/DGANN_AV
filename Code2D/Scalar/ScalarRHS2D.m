function [rhsQ] = ScalarRHS2D(Q,time,AdvectionVelocity)

% function [rhsQ] = EulerRHS2D(Q,time, ExactSolutionBC);
% Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
%             with a local Lax-Friedrich flux

Globals2D_DG;

vmapM = reshape(vmapM, Nfp*Nfaces, K); vmapP = reshape(vmapP, Nfp*Nfaces, K);

% 1. Compute volume contributions (NOW INDEPENDENT OF SURFACE TERMS)
gamma = 1.4;
[F,G] = ScalarFlux2D(Q, model, AdvectionVelocity);

% Compute weak derivatives
dFdr = Drw*F(:,:,1); dFds = Dsw*F(:,:,1);
dGdr = Drw*G(:,:,1); dGds = Dsw*G(:,:,1);
rhsQ(:,:,1) = (rx.*dFdr + sx.*dFds) + (ry.*dGdr + sy.*dGds);

    
% 2. Compute surface contributions 
% 2.1 evaluate '-' and '+' traces of conservative variables
QM(:,:,1) = Q(vmapM); QP(:,:,1) = Q(vmapP);


% 2.2 set boundary conditions by modifying positive traces
% if(~isempty(ExactSolutionBC))
%   QP = feval(ExactSolutionBC, Fx, Fy, nx, ny, mapI, mapO, mapW, mapC, QP, time);
% end

% 2.3 evaluate primitive variables & flux functions at '-' and '+' traces
[fM,gM] = ScalarFlux2D(QM, model, AdvectionVelocity);
[fP,gP] = ScalarFlux2D(QP, model, AdvectionVelocity);

% 2.4 Compute local Lax-Friedrichs/Rusonov numerical fluxes
lambda = max( Get_scalar_eig2D(QM,model,AdvectionVelocity),  ...
	          Get_scalar_eig2D(QP,model,AdvectionVelocity));
lambda = reshape(lambda, Nfp, Nfaces*K);
lambda = ones(Nfp, 1)*max(lambda, [], 1); 
lambda = reshape(lambda, Nfp*Nfaces, K);

% 2.5 Lift fluxes
flux = nx.*(fP(:,:,1) + fM(:,:,1)) + ny.*(gP(:,:,1) + gM(:,:,1)) + ...
       lambda.*(QM(:,:,1) - QP(:,:,1));
rhsQ(:,:,1) = rhsQ(:,:,1) - LIFT*(Fscale.*flux/2);

return;

