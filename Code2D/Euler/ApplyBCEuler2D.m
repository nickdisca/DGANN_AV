function QG = ApplyBCEuler2D(Q,time,Problem,Mesh)
 
% function  QG = ApplyBCEuler2D(Q,time,gas_gamma)
% Purpose: Set Boundary Condition in ghost elements. 
% NOTE: FOR INFLOW AND OUTFLOW TYPE BC, WE CURRENTLY MODIFY THE FACE
% VALUE ACCORDING TO SUBSONIC AND SUPERSONIC INFLOW/OUTFLOW 

rho  = Q(:,:,1);
rhou = Q(:,:,2);
rhov = Q(:,:,3);
Ener = Q(:,:,4);
pre  = Euler_Pressure2D(Q,Problem.gas_gamma);
u    = rhou./rho;
v    = rhov./rho;
a    = sqrt(Problem.gas_gamma*pre./rho);
rhoG = zeros(Mesh.Np,Mesh.KG);
uG   = zeros(Mesh.Np,Mesh.KG);
vG   = zeros(Mesh.Np,Mesh.KG);
preG = zeros(Mesh.Np,Mesh.KG);

BC_keys = Mesh.BC_ess_flags(:,1);

% Applying BC on triangles
for k=1:length(BC_keys)
    ckey   = BC_keys(k);
    bctype = Mesh.BC_ess_flags(k,2);
    
    [rhoBC,uBC,vBC,preBC] = BC(ckey,time,Problem.gas_gamma,Problem.gas_const);
    
    for f=1:3
        bc_ind = find(Mesh.BCTag(:,f) == ckey);
        gc_ind = Mesh.EToGE(bc_ind,f)';
        
        if(bctype==Mesh.BC_ENUM.Slip)
            nxa = mean(Mesh.nx(1+(f-1)*Mesh.Nfp:f*Mesh.Nfp,bc_ind'));
            nya = mean(Mesh.ny(1+(f-1)*Mesh.Nfp:f*Mesh.Nfp,bc_ind'));
            rhoG(:,gc_ind) = rho(Mesh.MMAP(:,f),bc_ind');
            uG(:,gc_ind)   = u(Mesh.MMAP(:,f),bc_ind') - 2*(u(Mesh.MMAP(:,f),bc_ind').*nxa + v(Mesh.MMAP(:,f),bc_ind').*nya).*nxa;
            vG(:,gc_ind)   = v(Mesh.MMAP(:,f),bc_ind') - 2*(u(Mesh.MMAP(:,f),bc_ind').*nxa + v(Mesh.MMAP(:,f),bc_ind').*nya).*nya;
            preG(:,gc_ind) = pre(Mesh.MMAP(:,f),bc_ind');
        elseif(bctype==Mesh.BC_ENUM.Dirichlet)
            rhoG(:,gc_ind) = rhoBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            uG(:,gc_ind)   = uBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            vG(:,gc_ind)   = vBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            preG(:,gc_ind) = preBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
%             rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
%             uG(:,gc_ind)   = u(MMAP(:,f),bc_ind');
%             vG(:,gc_ind)   = v(MMAP(:,f),bc_ind');
%             preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
        elseif( bctype==Mesh.BC_ENUM.In)
            rhoG(:,gc_ind) = rhoBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            uG(:,gc_ind)   = uBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            vG(:,gc_ind)   = vBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            preG(:,gc_ind) = preBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
%             rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
%             uG(:,gc_ind)   = u(MMAP(:,f),bc_ind');
%             vG(:,gc_ind)   = v(MMAP(:,f),bc_ind');
%             preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
        elseif(bctype == Mesh.BC_ENUM.Out)
%             rhoG(:,gc_ind) = rho(MMAP(:,f),bc_ind');
%             uG(:,gc_ind)   = u(MMAP(:,f),bc_ind');
%             vG(:,gc_ind)   = v(MMAP(:,f),bc_ind');
%             preG(:,gc_ind) = pre(MMAP(:,f),bc_ind');
            rhoG(:,gc_ind) = rhoBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            uG(:,gc_ind)   = uBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            vG(:,gc_ind)   = vBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind));
            preG(:,gc_ind) = preBC(Mesh.xG(:,gc_ind),Mesh.yG(:,gc_ind)); 
        else
            error('Unknown boundary condition for Euler equations');
        end 
    end
end

% Modifying traces for In and Out BC
rhoB = rhoG(Mesh.Fmask(:,Mesh.Nfaces),:);
uB   = uG(Mesh.Fmask(:,Mesh.Nfaces),:);
vB   = vG(Mesh.Fmask(:,Mesh.Nfaces),:);
preB = preG(Mesh.Fmask(:,Mesh.Nfaces),:);
rhoM = rho(Mesh.vmapM);
uM   = u(Mesh.vmapM);
vM   = v(Mesh.vmapM);
preM = pre(Mesh.vmapM);
for k=1:length(BC_keys)
    ckey   = BC_keys(k);
    bctype = Mesh.BC_ess_flags(k,2);
    mapBC  = Mesh.mapBC_list(ckey);
    %vmapBC = Mesh.vmapBC_list(ckey);
    gelem  = Mesh.GEBC_list(ckey);
    
    if(bctype==Mesh.BC_ENUM.In)
        [rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem)] =  EulerInflow2D(rhoM(mapBC),uM(mapBC),vM(mapBC),preM(mapBC),...
                                           rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem), Mesh.nx(mapBC), Mesh.ny(mapBC),Problem.gas_gamma);
    elseif(bctype==Mesh.BC_ENUM.Out)
        [rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem)] =  EulerOutflow2D(rhoM(mapBC),uM(mapBC),vM(mapBC),preM(mapBC),...
                                           rhoB(:,gelem),uB(:,gelem),vB(:,gelem),preB(:,gelem), Mesh.nx(mapBC), Mesh.ny(mapBC),Problem.gas_gamma);
    end
end

rhoG(Mesh.Fmask(:,Mesh.Nfaces),:) = rhoB;
uG(Mesh.Fmask(:,Mesh.Nfaces),:)   = uB;
vG(Mesh.Fmask(:,Mesh.Nfaces),:)   = vB;
preG(Mesh.Fmask(:,Mesh.Nfaces),:) = preB;
    

QG(:,:,1) = rhoG; QG(:,:,2) = rhoG.*uG; QG(:,:,3) = rhoG.*vG;
QG(:,:,4) = Euler_Energy2D(rhoG,uG,vG,preG,Problem.gas_gamma);

return;


