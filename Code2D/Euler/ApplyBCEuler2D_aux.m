function GG = ApplyBCEuler2D_aux(G,time,Problem,Mesh)
 
% function  QG = ApplyBCEuler2D(Q,time,gas_gamma)
% Purpose: Set Boundary Condition in ghost elements for gradient variable. 

G1  = G(:,:,1);
G2 = G(:,:,2);
G3 = G(:,:,3);
G4 = G(:,:,4);
G1G = zeros(Mesh.Np,Mesh.KG);
G2G   = zeros(Mesh.Np,Mesh.KG);
G3G   = zeros(Mesh.Np,Mesh.KG);
G4G = zeros(Mesh.Np,Mesh.KG);

BC_keys = Mesh.BC_ess_flags(:,1);

% Applying BC on triangles
for k=1:length(BC_keys)
    ckey   = BC_keys(k);
    bctype = Mesh.BC_ess_flags(k,2);
    
    for f=1:3
        bc_ind = find(Mesh.BCTag(:,f) == ckey);
        gc_ind = Mesh.EToGE(bc_ind,f)';
        
        if(bctype==Mesh.BC_ENUM.Slip)
            G1G(:,gc_ind) = -G1(MMAP(:,f),bc_ind');
            G2G(:,gc_ind)   = G2(MMAP(:,f),bc_ind');
            G3G(:,gc_ind)   = G3(MMAP(:,f),bc_ind');
            G4G(:,gc_ind) = -G4(MMAP(:,f),bc_ind');
        elseif(bctype==Mesh.BC_ENUM.Dirichlet)
            G1G(:,gc_ind) = G1(MMAP(:,f),bc_ind');
            G2G(:,gc_ind)   = G2(MMAP(:,f),bc_ind');
            G3G(:,gc_ind)   = G3(MMAP(:,f),bc_ind');
            G4G(:,gc_ind) = G4(MMAP(:,f),bc_ind');
        elseif( bctype==Mesh.BC_ENUM.In)
            G1G(:,gc_ind) = G1(MMAP(:,f),bc_ind');
            G2G(:,gc_ind)   = G2(MMAP(:,f),bc_ind');
            G3G(:,gc_ind)   = G3(MMAP(:,f),bc_ind');
            G4G(:,gc_ind) = G4(MMAP(:,f),bc_ind');
        elseif(bctype == Mesh.BC_ENUM.Out)
            G1G(:,gc_ind) = -G1(MMAP(:,f),bc_ind');
            G2G(:,gc_ind)   = -G2(MMAP(:,f),bc_ind');
            G3G(:,gc_ind)   = -G3(MMAP(:,f),bc_ind');
            G4G(:,gc_ind) = -G4(MMAP(:,f),bc_ind'); 
        else
            error('Unknown boundary condition for Euler equations');
        end 
    end
end


GG(:,:,1) = G1G; GG(:,:,2) = G2G; GG(:,:,3) = G3G;
GG(:,:,4) = G4G;


return;
