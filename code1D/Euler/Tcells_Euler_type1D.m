function ind = Tcells_Euler_type1D(q,Problem,Mesh,Limit,Net)

% Purpose: find all the troubled-cells for variable(s)


q_ext(:,:,1)  = Apply_BC1D(q(:,:,1),Problem.bc_cond(1,:));
q_ext(:,:,2)  = Apply_BC1D(q(:,:,2),Problem.bc_cond(2,:));
q_ext(:,:,3)  = Apply_BC1D(q(:,:,3),Problem.bc_cond(3,:));

switch Limit.ind_var
    case 'density'
        ind = Find_Tcells1D(q_ext(:,:,1),Mesh,Limit,Net);
    case 'velocity'
        ind = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1),Mesh,Limit,Net);
    case 'pressure'
        pre = (Problem.gas_gamma-1)*(q_ext(:,:,3) - 0.5*q_ext(:,:,2).^2./q_ext(:,:,1)); 
        ind = Find_Tcells1D(pre,Mesh,Limit,Net);  
    case 'prim'
        pre = (Problem.gas_gamma-1)*(q_ext(:,:,3) - 0.5*q_ext(:,:,2).^2./q_ext(:,:,1)); 
        ind_d = Find_Tcells1D(q_ext(:,:,1),Mesh,Limit,Net);
        ind_v = Find_Tcells1D(q_ext(:,:,2)./q_ext(:,:,1),Mesh,Limit,Net);
        ind_p = Find_Tcells1D(pre,Mesh,Limit,Net);
        ind   = unique([ind_d,ind_v,ind_p]);
    case 'con'
        ind_d = Find_Tcells1D(q_ext(:,:,1),Mesh,Limit,Net);
        ind_m = Find_Tcells1D(q_ext(:,:,2),Mesh,Limit,Net);
        ind_e = Find_Tcells1D(q_ext(:,:,3),Mesh,Limit,Net);
        ind   = unique([ind_d,ind_m,ind_e]);    
    otherwise
        error('Unknown indicator variable %s',Limit.ind_var)    
end



return
