function ind = Find_Tcells2D(Q,QG,Limit,Mesh,Net)

if(strcmp(Limit.Indicator,'NONE'))
   ind = [];
   return   
elseif(strcmp(Limit.Indicator,'ALL'))
   ind = 1:Mesh.K;
   return      
end

if(Limit.Filter_const)
    nconst_ind = FindNonConstCells2D(Q);
end

switch Limit.Indicator
    
    case 'TVB'
        ind = TVB_Indicator2D(Q,QG,Mesh,Limit.TVBM,Limit.TVBnu); 
        ind = intersect(ind,nconst_ind); % Need to make this a pre-processing step
        
    case 'NN'
        ind = NN_Indicator2D(Q,QG,Mesh,Net,nconst_ind);  
        
    otherwise
        error('Unknown indicator type %s',Limit.Indicator)
end




% if(Limit.Remove_iso)
% %     niso_ind = FindNonIsoFlaggedCells2D(ind);
% % else
% %     nconst_ind = 1:Mesh.K;
%    error('To DO');
% end

        
return
