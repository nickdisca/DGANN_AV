function ind = Find_Tcells2D(Q,Indicator,Limiter)

if(strcmp(Limiter,'NONE'))
   ind = [];
   return
end
switch Indicator
    
    case 'NONE'
        ind = 1:length(Q(1,:));
    
    case 'TVB'
        ind = TVB_Indicator2D(Q);
        
    case 'TVB2'
        ind = TVB_Indicator2D_ncon(Q);    
        
    case 'NN'
        ind = NN_Indicator2D(Q);
        
    case 'NN_Pwise'
        ind = NN_Indicator2D(Q);
        
    case 'NN_modal_Pwise'
        ind = NN_Indicator2D_modal(Q);    
        
    case 'NN_modal_patch_Pwise'
        ind = NN_Indicator2D_modal_patch(Q);     
        
    otherwise
        error('Unknown indicator type %s',Indicator)
end
        
return
