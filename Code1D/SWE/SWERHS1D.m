function [rhsdepth, rhsdvel] = SWERHS1D(depth, dvel, bc_type, gravity)

% function [rhsdepth, rhsdvel] = SWERHS1D(depth, dvel, bc_type, gravity)
% Purpose  : Evaluate RHS flux in 1D Shallow water

Globals1D;

% compute maximum velocity for LF flux
cvel = sqrt(gravity*depth); lm = abs(dvel./depth)+cvel;

% Compute fluxes
depthf = dvel; dvelf= 0.5*gravity*depth.^2 + dvel.^2./depth;

% Compute jumps at internal faces
ddepth  =zeros(Nfp*Nfaces,K);  ddepth(:)  = depth(vmapM) -  depth(vmapP); 
ddvel   =zeros(Nfp*Nfaces,K);  ddvel(:)   = dvel(vmapM)  - dvel(vmapP);
ddepthf =zeros(Nfp*Nfaces,K);  ddepthf(:) = depthf(vmapM)-  depthf(vmapP); 
ddvelf  =zeros(Nfp*Nfaces,K);  ddvelf(:)  = dvelf(vmapM) - dvelf(vmapP);
LFc     =zeros(Nfp*Nfaces,K);  LFc(:)     = max(lm(vmapP),lm(vmapM));


% Compute fluxes at interfaces
ddepthf(:) = nx(:).*ddepthf(:)/2.0-LFc(:)/2.0.*ddepth(:); 
ddvelf(:)  = nx(:).*ddvelf(:)/2.0-LFc(:)/2.0.*ddvel(:); 


% Boundary conditions
[depthin,dvelin,depthout,dvelout] = SWEBC1D(bc_type,depth(vmapI),depth(vmapO),...
                                             dvel(vmapI),dvel(vmapO));

% Set fluxes at inflow/outflow
depthfin =dvelin; dvelfin=0.5*gravity*depthin.^2 + dvelin.^2./depthin;
lmI=lm(vmapI)/2; nxI=nx(mapI);
ddepthf (mapI)=nxI*(depthf (vmapI)-depthfin )/2.0-lmI*(depth(vmapI) -depthin);  
ddvelf(mapI)  =nxI*(dvelf(vmapI)-dvelfin)/2.0-lmI*(dvel(vmapI)-dvelin);

depthfout =dvelout; dvelfout=0.5*gravity*depthout.^2 + dvelout.^2./depthout;
lmO=lm(vmapO)/2; nxO=nx(mapO);
ddepthf (mapO)=nxO*(depthf (vmapO)-depthfout )/2.0-lmO*(depth(vmapO) -depthout);  
ddvelf(mapO)  =nxO*(dvelf(vmapO)-dvelfout)/2.0-lmO*(dvel(vmapO)-dvelout);

% compute right hand sides of the PDE's
rhsdepth  = -rx.*(Dr*depthf)  + LIFT*(Fscale.*ddepthf);
rhsdvel   = -rx.*(Dr*dvelf) + LIFT*(Fscale.*ddvelf);
return
