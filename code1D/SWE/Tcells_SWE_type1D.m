function ind = Tcells_SWE_type1D(q,ind_var,bc_cond)

% Purpose: find all the troubled-cells for variable(s)

q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),bc_cond(2,:));

switch ind_var
    case 'depth'
        ind = Find_Tcells1D(q_ext(:,:,1));
    case 'velocity'
        ind = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1));
    case 'dv'
        ind_d = Find_Tcells1D(q_ext(:,:,1));
        ind_v = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1));
        ind   = unique([ind_d,ind_v]);
    otherwise
        error('Unknown indicator variable %s',ind_var)    
end



return
