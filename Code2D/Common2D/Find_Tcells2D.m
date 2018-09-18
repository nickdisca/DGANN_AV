function ind = Find_Tcells2D(Q,QG,Indicator,Limiter)

if(strcmp(Limiter,'NONE'))
   ind = [];
   return
end
switch Indicator
    
    case 'NONE'
        ind = 1:length(Q(1,:));
    
    case 'TVB'
        ind = TVB_Indicator2D(Q,QG);
        
    case 'TVB2'
        ind = TVB_Indicator2D_ncon(Q,QG);    
        
    case 'NN'
        ind = NN_Indicator2D(Q,QG);
        
    case 'NN_Pwise'
        ind = NN_Indicator2D(Q,QG);
        
    case 'NN_modal_Pwise'
        ind = NN_Indicator2D_modal(Q,QG);    
        
    case 'NN_modal_patch_Pwise'
        ind = NN_Indicator2D_modal_patch(Q,QG);     
        
    otherwise
        error('Unknown indicator type %s',Indicator)
end
        
return
