function [rhsQ] = EulerRHS2D(Q,time,mu,Problem,Mesh)

% function [rhsQ] = EulerRHS2D(Q,time, gas_gamma, gas_const);
% Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
%             with a local Lax-Friedrich flux

% Apply boundary conditions
QG = ApplyBCEuler2D(Q,time,Problem,Mesh);

% Solution traces
QM = zeros(numel(Mesh.Fmask),Mesh.K,4);
QP = zeros(numel(Mesh.Fmask),Mesh.K,4);
for n=1:4
  Qn = Q(:,:,n);
  QM(:,:,n) = Qn(Mesh.vmapM); QP(:,:,n) = Qn(Mesh.vmapP);
end

% Set boundary conditions by modifying positive traces
if(~isempty(Mesh.mapBC_list))
    for n=1:4
        QPn = QP(:,:,n);
        QPn(Mesh.mapB) = QG(Mesh.Fmask(:,Mesh.Nfaces),:,n);
        QP(:,:,n) = QPn;
    end
end

% Compute numerical fluxes of auxiliary equation
bd_flux_X=zeros(size(Q)); bd_flux_Y=zeros(size(Q));
for n=1:4
    bd_flux_X(:,:,n)=Mesh.LIFT*(Mesh.Fscale.*((QP(:,:,n)+QM(:,:,n))/2.*Mesh.nx));
    bd_flux_Y(:,:,n)=Mesh.LIFT*(Mesh.Fscale.*((QP(:,:,n)+QM(:,:,n))/2.*Mesh.ny));
end

%Compute internal fluxes for auxiliary variable
int_flux_X=zeros(size(Q));
int_flux_Y=zeros(size(Q));
for n=1:4
    int_flux_X(:,:,n)=Mesh.rx.*(Mesh.Drw*Q(:,:,n))+Mesh.sx.*(Mesh.Dsw*Q(:,:,n));
    int_flux_Y(:,:,n)=Mesh.ry.*(Mesh.Drw*Q(:,:,n))+Mesh.sy.*(Mesh.Dsw*Q(:,:,n));
end

% Compute auxiliary variables with viscosity
gX=zeros(size(Q));
gY=zeros(size(Q));
for n=1:4
    gX(:,:,n)=mu(:,:,n).*(bd_flux_X(:,:,n)-int_flux_X(:,:,n));
    gY(:,:,n)=mu(:,:,n).*(bd_flux_Y(:,:,n)-int_flux_Y(:,:,n));
end

% Apply boundary conditions to auxiliary variable
gXG = ApplyBCEuler2D_aux(gX,time,Problem,Mesh);
gYG = ApplyBCEuler2D_aux(gY,time,Problem,Mesh);

% Variable traces
gXM = zeros(numel(Mesh.Fmask),Mesh.K,4);
gXP = zeros(numel(Mesh.Fmask),Mesh.K,4);
gYM = zeros(numel(Mesh.Fmask),Mesh.K,4);
gYP = zeros(numel(Mesh.Fmask),Mesh.K,4);
for n=1:4
  Gn = gX(:,:,n); gXM(:,:,n) = Gn(Mesh.vmapM); gXP(:,:,n) = Gn(Mesh.vmapP);
  Gn = gY(:,:,n); gYM(:,:,n) = Gn(Mesh.vmapM); gYP(:,:,n) = Gn(Mesh.vmapP);
end

% Set boundary conditions by modifying positive traces
if(~isempty(Mesh.mapBC_list))
    for n=1:4
        gXPn = gXP(:,:,n);
        gXPn(Mesh.mapB) = gXG(Mesh.Fmask(:,Mesh.Nfaces),:,1);
        gXP(:,:,n) = gXPn;
        
        gYPn = gYP(:,:,n);
        gYPn(Mesh.mapB) = gYG(Mesh.Fmask(:,Mesh.Nfaces),:,1);
        gYP(:,:,n) = gYPn;
    end
end



% Compute volume contributions 
[F,G,~,~,~,~] = EulerExactFluxes2D(Q, Problem.gas_gamma);

% Compute weak derivatives
dFdr=zeros(size(Q));
dFds=zeros(size(Q));
dGdr=zeros(size(Q));
dGds=zeros(size(Q));
for n=1:4
  dFdr(:,:,n) = Mesh.Drw*F(:,:,n); dFds(:,:,n) = Mesh.Dsw*F(:,:,n);
  dGdr(:,:,n) = Mesh.Drw*G(:,:,n); dGds(:,:,n) = Mesh.Dsw*G(:,:,n);
end

int_flux=zeros(size(Q));
for n=1:4
    % Compute volume term of physical variable
    int_flux(:,:,n) = (Mesh.rx.*dFdr(:,:,n) + Mesh.sx.*dFds(:,:,n)) + (Mesh.ry.*dGdr(:,:,n) + Mesh.sy.*dGds(:,:,n));
    % Volume term for auxiliary variable
    int_flux(:,:,n)=int_flux(:,:,n)-(Mesh.rx.*(Mesh.Drw*gX(:,:,n))+Mesh.sx.*(Mesh.Dsw*gX(:,:,n))+Mesh.ry.*(Mesh.Drw*gY(:,:,n))+Mesh.sy.*(Mesh.Dsw*gY(:,:,n)));
end
    


% 2.3 evaluate primitive variables & flux functions 
[fM,gM,rhoM,uM,vM,pM] = EulerExactFluxes2D(QM, Problem.gas_gamma);
[fP,gP,rhoP,uP,vP,pP] = EulerExactFluxes2D(QP, Problem.gas_gamma);

% 2.4 Compute maximum wave speed on edges
lambda = max( sqrt(uM.^2+vM.^2) + sqrt(abs(Problem.gas_gamma*pM./rhoM)),  ...
	      sqrt(uP.^2+vP.^2) + sqrt(abs(Problem.gas_gamma*pP./rhoP)));
lambda = reshape(lambda, Mesh.Nfp, Mesh.Nfaces*Mesh.K);
lambda = ones(Mesh.Nfp, 1)*max(lambda, [], 1); 
lambda = reshape(lambda, Mesh.Nfp*Mesh.Nfaces, Mesh.K);


bd_flux=zeros(size(Q));
for n=1:4
    % Boundary flux for physical variable (Lax-Friedrichs)
    nflux = Mesh.nx.*(fP(:,:,n) + fM(:,:,n)) + Mesh.ny.*(gP(:,:,n) + gM(:,:,n)) + ...
        lambda.*(QM(:,:,n) - QP(:,:,n));
    % Boundary flux from auxiliary variable
    nflux = nflux - ((gXP(:,:,n)+gXM(:,:,n)).*Mesh.nx+(gYP(:,:,n)+gYM(:,:,n)).*Mesh.ny);
    % Lift
    bd_flux(:,:,n) = Mesh.LIFT*(Mesh.Fscale.*nflux/2);
end

% Compute rhs   
rhsQ = int_flux - bd_flux;


return;




