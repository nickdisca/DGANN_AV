function ind = Tcells_Euler_type1D(q,gas_gamma, gas_const,ind_var,bc_cond)

% Purpose: find all the troubled-cells for variable(s)


q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),bc_cond(2,:));
q_ext(:,:,3)  = Apply_BC1D(q(:,:,3),bc_cond(3,:));

switch ind_var
    case 'density'
        ind = Find_Tcells1D(q_ext(:,:,1));
    case 'velocity'
        ind = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1));
    case 'pressure'
        pre = (gas_gamma-1)*(q_ext(:,:,3) - 0.5*q_ext(:,:,2).^2./q_ext(:,:,1)); 
        ind = Find_Tcells1D(pre);    
    case 'prim'
        pre = (gas_gamma-1)*(q_ext(:,:,3) - 0.5*q_ext(:,:,2).^2./q_ext(:,:,1)); 
        ind_d = Find_Tcells1D(q_ext(:,:,1));
        ind_v = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1));
        ind_p = Find_Tcells1D(pre);  
        ind   = unique([ind_d,ind_v,ind_p]);
    case 'energy' 
        ind = Find_Tcells1D(q_ext(:,:,3));
    case 'de' 
        ind_d = Find_Tcells1D(q_ext(:,:,1)); 
        ind_e = Find_Tcells1D(q_ext(:,:,3)); 
        ind   = unique([ind_d,ind_e]);
    otherwise
        error('Unknown indicator variable %s',ind_var)    
end



return
