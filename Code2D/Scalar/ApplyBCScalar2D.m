function QG = ApplyBCScalar2D(Q,time)
 
% function  QG = ApplyBCScalar2D(Q,time,gas_gamma)
% Purpose: Set Boundary Condition in ghost elements. 

Globals2D_DG;

QG = zeros(Np,KG,1);

BC_keys = BC_ess_flags(:,1);

% Applying BC on triangles
for k=1:length(BC_keys)
    ckey   = BC_keys(k);
    bctype = BC_ess_flags(k,2);
    
    QBC = BC(ckey,time);
    
    for f=1:3
        bc_ind = find(BCTag(:,f) == ckey);
        gc_ind = EToGE(bc_ind,f)';
        
        if(bctype==Sym)
            QG(:,gc_ind,1) = Q(MMAP(:,f),bc_ind',1);
        elseif(bctype==Dirichlet)
            QG(:,gc_ind,1) = QBC(xG(:,gc_ind),yG(:,gc_ind),1);  
        end 
    end
end


return;


