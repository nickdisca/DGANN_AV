function ind = Tcells_SWE_type1D(q,bc_cond,Mesh,Limit,Net)

% Purpose: find all the troubled-cells for variable(s)

q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),bc_cond(2,:));

switch Limit.ind_var
    
    case 'prim'
        ind_d = Find_Tcells1D(q_ext(:,:,1),Mesh,Limit,Net);
        ind_v = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1),Mesh,Limit,Net);
        ind   = unique([ind_d,ind_v]);
    
    case 'con'
        ind_d   = Find_Tcells1D(q_ext(:,:,1),Mesh,Limit,Net);
        ind_dis = Find_Tcells1D(q_ext(:,:,2),Mesh,Limit,Net);
        ind     = unique([ind_d,ind_dis]);    
        
    otherwise
        error('Unknown indicator variable %s',Limit.ind_var)    
end



return
