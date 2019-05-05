function ind = Tcells_Euler_type2D(Q,QG,gas_gamma,Limit,Mesh,Net)

% Purpose: find all the troubled-cells for variable(s)

if(strcmp(Limit.Indicator,'NONE'))
   ind = [];
   return   
elseif(strcmp(Limit.Indicator,'ALL'))
   ind = 1:Mesh.K;
   return      
end


switch Limit.ind_var
    case 'density'
        ind = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
    case 'velocity'
        ind_u = Find_Tcells2D(Q(:,:,2)./Q(:,:,1),QG(:,:,2)./QG(:,:,1),Limit,Mesh,Net);
        ind_v = Find_Tcells2D(Q(:,:,3)./Q(:,:,1),QG(:,:,3)./QG(:,:,1),Limit,Mesh,Net);
        ind   = unique([ind_u,ind_v]);
    case 'pressure'
        pre  = Euler_Pressure2D(Q,gas_gamma); 
        preG = Euler_Pressure2D(QG,gas_gamma); 
        ind  = Find_Tcells2D(pre,preG,Limit,Mesh,Net);    
    case 'prim'
        pre   = Euler_Pressure2D(Q,gas_gamma); 
        preG  = Euler_Pressure2D(QG,gas_gamma); 
        ind_d = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Indicator,Limiter);
        ind_u = Find_Tcells2D(Q(:,:,2)./Q(:,:,1),QG(:,:,2)./QG(:,:,1),Indicator,Limiter);
        ind_v = Find_Tcells2D(Q(:,:,3)./Q(:,:,1),QG(:,:,3)./QG(:,:,1),Indicator,Limiter);
        ind_p = Find_Tcells2D(pre,preG,Limit,Mesh,Net);  
        ind   = unique([ind_d,ind_u,ind_v,ind_p]);
    case 'con'
        ind_1 = Find_Tcells2D(Q(:,:,1),QG(:,:,1),Limit,Mesh,Net);
        ind_2 = Find_Tcells2D(Q(:,:,2),QG(:,:,2),Limit,Mesh,Net);
        ind_3 = Find_Tcells2D(Q(:,:,3),QG(:,:,3),Limit,Mesh,Net);
        ind_4 = Find_Tcells2D(Q(:,:,4),QG(:,:,4),Limit,Mesh,Net);  
        ind   = unique([ind_1,ind_2,ind_3,ind_4]);    
    otherwise
        error('Unknown indicator variable %s',ind_var)    
end

return
