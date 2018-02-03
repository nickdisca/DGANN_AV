function [rhsu] = ScalarRHS1D_strong(u,flux,dflux)

% TO DO-----------


% Purpose  : Evaluate RHS flux in 1D Scalar equation using the strong form

Globals1D_DG;

lm = abs(dflux(u));
% Define field differences at faces
du  = zeros(Nfp*Nfaces,K); du(:)  = u(vmapM)-u(vmapP);
duf = zeros(Nfp*Nfaces,K); duf(:) = flux(u(vmapM))-flux(u(vmapP));
LFc = zeros(Nfp*Nfaces,K); LFc(:) = max(lm(vmapP),lm(vmapM));

% Compute fluxes at interfaces
duf(:) =nx(:).*duf(:)/2.0-LFc(:)/2.0.*du(:); 

% impose boundary condition 
[uin,uout] = Scalar_BC1D(bc_type,u(vmapI),u(vmapO));
lmI=lm(vmapI)/2; nxI=nx(mapI);
lmO=lm(vmapO)/2; nxO=nx(mapO);
duf(mapI)=nxI*(flux(u(vmapI))-flux(uin))/2.0-lmI*(u(vmapI) -uin); 
duf(mapO)=nxO*(flux(u(vmapO))-flux(uout))/2.0-lmO*(u(vmapO) -uout); 

% compute right hand sides of the semi-discrete PDE
rhsu  = -rx.*(Dr*flux(u))  + LIFT*(Fscale.*duf);
return
