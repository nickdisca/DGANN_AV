function QG = ApplyBCEuler2D(Q,time,gas_gamma)
 
% function  QG = ApplyBCEuler2D(Q,time,gas_gamma)
% Purpose: Set Boundary Condition in ghost elements. 
% NOTE: FOR INFLOW AND OUTFLOW TYPE BC, WE CURRENTLY MODIFY THE FACE
% VALUE ACCORDING TO SUBSONIC AND SUPERSONIC INFLOW/OUTFLOW 

Globals2D_DG;

rho  = Q(:,:,1);
rhou = Q(:,:,2);
rhov = Q(:,:,3);
Ener = Q(:,:,4);
pre  = Euler_Pressure2D(Q,gas_gamma);
u    = rhou./rho;
v    = rhov./rho;
a    = sqrt(gas_gamma*pre./rho);
rhoG = zeros(Np,KG);
uG   = zeros(Np,KG);
vG   = zeros(Np,KG);
preG = zeros(Np,KG);

BC_keys = BC_ess_flags(:,1);

% Applying BC on triangles
for k=1:length(BC_keys)
    ckey   = BC_keys(k);
    bctype = BC_ess_flags(k,2);
    
    [rhoBC,uBC,vBC,preBC] = BC(ckey,time,gas_gamma);
    
    for f=1:3
        bc_ind = find(BCTag(:,f) == ckey);
        gc_ind = EToGE(bc_ind,f)';
        
        if(bctype==Slip)
            nxa = mean(nx(1+(f-1)*Nfp:f*Nfp,bc_ind'));
            nya = mean(ny(1+(f-1)*Nfp:f*Nfp,bc_ind'));
            rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
            uG(:,gc_ind)   = u(MMAP(:,f),bc_ind') - 2*(u(MMAP(:,f),bc_ind').*nxa + v(MMAP(:,f),bc_ind').*nya).*nxa;
            vG(:,gc_ind)   = v(MMAP(:,f),bc_ind') - 2*(u(MMAP(:,f),bc_ind').*nxa + v(MMAP(:,f),bc_ind').*nya).*nya;
            preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
        elseif(bctype==Dirichlet)
            rhoG(:,gc_ind) = rhoBC(xG(:,gc_ind),yG(:,gc_ind));
            uG(:,gc_ind)   = uBC(xG(:,gc_ind),yG(:,gc_ind));
            vG(:,gc_ind)   = vBC(xG(:,gc_ind),yG(:,gc_ind));
            preG(:,gc_ind) = preBC(xG(:,gc_ind),yG(:,gc_ind));
%             rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
%             uG(:,gc_ind)   = u(MMAP(:,f),bc_ind');
%             vG(:,gc_ind)   = v(MMAP(:,f),bc_ind');
%             preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
        elseif( bctype==In)
            rhoG(:,gc_ind) = rhoBC(xG(:,gc_ind),yG(:,gc_ind));
            uG(:,gc_ind)   = uBC(xG(:,gc_ind),yG(:,gc_ind));
            vG(:,gc_ind)   = vBC(xG(:,gc_ind),yG(:,gc_ind));
            preG(:,gc_ind) = preBC(xG(:,gc_ind),yG(:,gc_ind));
%             rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
%             uG(:,gc_ind)   = u(MMAP(:,f),bc_ind');
%             vG(:,gc_ind)   = v(MMAP(:,f),bc_ind');
%             preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
        elseif(bctype == Out)
%             rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
%             uG(:,gc_ind)   = u(MMAP(:,f),bc_ind');
%             vG(:,gc_ind)   = v(MMAP(:,f),bc_ind');
%             preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
            rhoG(:,gc_ind) = rhoBC(xG(:,gc_ind),yG(:,gc_ind));
            uG(:,gc_ind)   = uBC(xG(:,gc_ind),yG(:,gc_ind));
            vG(:,gc_ind)   = vBC(xG(:,gc_ind),yG(:,gc_ind));
            preG(:,gc_ind) = preBC(xG(:,gc_ind),yG(:,gc_ind));    
        end 
    end
end

% Modifying traces for In and Out BC
rhoB = rhoG(Fmask(:,Nfaces),:);
uB   = uG(Fmask(:,Nfaces),:);
vB   = vG(Fmask(:,Nfaces),:);
preB = preG(Fmask(:,Nfaces),:);
rhoM = rho(vmapM);
uM   = u(vmapM);
vM   = v(vmapM);
preM = pre(vmapM);
for k=1:length(BC_keys)
    ckey   = BC_keys(k);
    bctype = BC_ess_flags(k,2);
    mapBC  = mapBC_list(ckey);
    vmapBC = vmapBC_list(ckey);
    gelem  = GEBC_list(ckey);
    
    if(bctype==In)
        [rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem)] =  EulerInflow2D(rhoM(mapBC),uM(mapBC),vM(mapBC),preM(mapBC),...
                                           rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem), nx(mapBC), ny(mapBC),gas_gamma);
    elseif(bctype==Out)
        [rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem)] =  EulerOutflow2D(rhoM(mapBC),uM(mapBC),vM(mapBC),preM(mapBC),...
                                           rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem), nx(mapBC), ny(mapBC),gas_gamma);
    end
end

rhoG(Fmask(:,Nfaces),:) = rhoB;
uG(Fmask(:,Nfaces),:)   = uB;
vG(Fmask(:,Nfaces),:)   = vB;
preG(Fmask(:,Nfaces),:) = preB;
    

QG(:,:,1) = rhoG; QG(:,:,2) = rhoG.*uG; QG(:,:,3) = rhoG.*vG;
QG(:,:,4) = Euler_Energy2D(rhoG,uG,vG,preG,gas_gamma);

return;


