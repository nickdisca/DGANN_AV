function ind = Tcells_SWE_type(q,ind_var)

% Purpose: find all the troubled-cells for variable(s)

switch ind_var
    case 'depth'
        ind = Find_Tcells(q(:,:,1));
    case 'velocity'
        ind = Find_Tcells(q(:,:,2)./q(:,:,1));
    case 'dv'
        ind_d = Find_Tcells(q(:,:,1));
        ind_v = Find_Tcells(q(:,:,2)./q(:,:,1));
        ind   = unique([ind_d,ind_v]);
    otherwise
        error('Unknown indicator variable %s',ind_var)    
end



return
