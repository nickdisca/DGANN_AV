function [rhsQ] = EulerRHS2D(Q,time,Problem,Mesh)

% function [rhsQ] = EulerRHS2D(Q,time, gas_gamma, gas_const);
% Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
%             with a local Lax-Friedrich flux

% Apply boundary conditions
QG = ApplyBCEuler2D(Q,time,Problem,Mesh);

% 1. Compute volume contributions (NOW INDEPENDENT OF SURFACE TERMS)
[F,G,~,~,~,~] = EulerExactFluxes2D(Q, Problem.gas_gamma);

% Compute weak derivatives
rhsQ = zeros(size(Q));
for n=1:4
  dFdr = Mesh.Drw*F(:,:,n); dFds = Mesh.Dsw*F(:,:,n);
  dGdr = Mesh.Drw*G(:,:,n); dGds = Mesh.Dsw*G(:,:,n);
  rhsQ(:,:,n) = (Mesh.rx.*dFdr + Mesh.sx.*dFds) + (Mesh.ry.*dGdr + Mesh.sy.*dGds);
end
    
% 2. Compute surface contributions 
% 2.1 evaluate '-' and '+' traces of conservative variables
QM = zeros(numel(Mesh.Fmask),Mesh.K,4);
QP = zeros(numel(Mesh.Fmask),Mesh.K,4);
for n=1:4
  Qn = Q(:,:,n);
  QM(:,:,n) = Qn(Mesh.vmapM); QP(:,:,n) = Qn(Mesh.vmapP);
end

% % 2.2 set boundary conditions by modifying positive traces
if(~isempty(Mesh.mapBC_list))
    for n=1:4
        QPn = QP(:,:,n);
        QPn(Mesh.mapB) = QG(Mesh.Fmask(:,Mesh.Nfaces),:,n);
        QP(:,:,n) = QPn;
    end
end

% 2.3 evaluate primitive variables & flux functions at '-' and '+' traces
[fM,gM,rhoM,uM,vM,pM] = EulerExactFluxes2D(QM, Problem.gas_gamma);
[fP,gP,rhoP,uP,vP,pP] = EulerExactFluxes2D(QP, Problem.gas_gamma);

% 2.4 Compute local Lax-Friedrichs/Rusonov numerical fluxes
lambda = max( sqrt(uM.^2+vM.^2) + sqrt(abs(Problem.gas_gamma*pM./rhoM)),  ...
	      sqrt(uP.^2+vP.^2) + sqrt(abs(Problem.gas_gamma*pP./rhoP)));
lambda = reshape(lambda, Mesh.Nfp, Mesh.Nfaces*Mesh.K);
lambda = ones(Mesh.Nfp, 1)*max(lambda, [], 1); 
lambda = reshape(lambda, Mesh.Nfp*Mesh.Nfaces, Mesh.K);

% 2.5 Lift fluxes
for n=1:4
  nflux = Mesh.nx.*(fP(:,:,n) + fM(:,:,n)) + Mesh.ny.*(gP(:,:,n) + gM(:,:,n)) + ...
      lambda.*(QM(:,:,n) - QP(:,:,n));
  rhsQ(:,:,n) = rhsQ(:,:,n) - Mesh.LIFT*(Mesh.Fscale.*nflux/2);
end
return;

