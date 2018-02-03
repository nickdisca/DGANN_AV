function ind = Tcells_Euler_type(q,gas_gamma, gas_const,ind_var)

% Purpose: find all the troubled-cells for variable(s)

switch ind_var
    case 'density'
        ind = Find_Tcells(q(:,:,1));
    case 'velocity'
        ind = Find_Tcells(q(:,:,2)./q(:,:,1));
    case 'pressure'
        pre = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1)); 
        ind = Find_Tcells(pre);    
    case 'prim'
        pre = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1)); 
        ind_d = Find_Tcells(q(:,:,1));
        ind_v = Find_Tcells(q(:,:,2)./q(:,:,1));
        ind_p = Find_Tcells(pre);  
        ind   = unique([ind_d,ind_v,ind_p]);
    case 'energy' 
        ind = Find_Tcells(q(:,:,3));
    case 'de' 
        ind_d = Find_Tcells(q(:,:,1)); 
        ind_e = Find_Tcells(q(:,:,3)); 
        ind   = unique([ind_d,ind_e]);
    otherwise
        error('Unknown indicator variable %s',ind_var)    
end



return
