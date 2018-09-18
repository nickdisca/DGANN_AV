function ind = Tcells_Euler_type2D(Q,QG,gas_gamma, gas_const,Indicator,ind_var,Limiter)

% Purpose: find all the troubled-cells for variable(s)

switch ind_var
    case 'density'
        ind = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter);
    case 'velocity'
        ind_u = Find_Tcells2D(Q(:,:,2)./Q(:,:,1),QG(:,:,2)./QG(:,:,1),Indicator,Limiter);
        ind_v = Find_Tcells2D(Q(:,:,3)./Q(:,:,1),QG(:,:,3)./QG(:,:,1),Indicator,Limiter);
        ind   = unique([ind_u,ind_v]);
    case 'pressure'
        pre  = Euler_Pressure2D(Q,gas_gamma); 
        preG = Euler_Pressure2D(QG,gas_gamma); 
        ind  = Find_Tcells2D(pre,preG,Indicator,Limiter);    
    case 'prim'
        pre   = Euler_Pressure2D(Q,gas_gamma); 
        preG  = Euler_Pressure2D(QG,gas_gamma); 
        ind_d = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter);
        ind_u = Find_Tcells2D(Q(:,:,2)./Q(:,:,1),QG(:,:,2)./QG(:,:,1),Indicator,Limiter);
        ind_v = Find_Tcells2D(Q(:,:,3)./Q(:,:,1),QG(:,:,3)./QG(:,:,1),Indicator,Limiter);
        ind_p = Find_Tcells2D(pre,preG,Indicator,Limiter);  
        ind   = unique([ind_d,ind_u,ind_v,ind_p]);
    case 'con'
        ind_1 = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter);
        ind_2 = Find_Tcells2D(Q(:,:,2),QG(:,:,2),Indicator,Limiter);
        ind_3 = Find_Tcells2D(Q(:,:,3),QG(:,:,3),Indicator,Limiter);
        ind_4 = Find_Tcells2D(Q(:,:,4),QG(:,:,4),Indicator,Limiter);  
        ind   = unique([ind_1,ind_2,ind_3,ind_4]);    
    case 'energy' 
        ind = Find_Tcells2D(Q(:,:,4),QG(:,:,4),Indicator,Limiter);
    case 'de' 
        ind_d = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter); 
        ind_e = Find_Tcells2D(Q(:,:,4),QG(:,:,4),Indicator,Limiter); 
        ind   = unique([ind_d,ind_e]);
    otherwise
        error('Unknown indicator variable %s',ind_var)    
end

return
