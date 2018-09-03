% Create save file base name (saved in th same folder)

data_fname = sprintf('%s2D_%s_P%d_%s',model,test_name,N,mtail);

if(strcmp(Limiter,'NONE'))
    data_fname = sprintf('%s_NONE',data_fname);
else
    if(strcmp(Indicator,'NONE'))
        data_fname = sprintf('%s_IND_%s',data_fname,Indicator);
    elseif(strcmp(Indicator,'TVB'))
        data_fname = sprintf('%s_IND_%s_%d',data_fname,Indicator,TVBM);
    elseif(strcmp(Indicator,'TVB2'))
        data_fname = sprintf('%s_IND_%s_%d',data_fname,Indicator,TVBM);     
    elseif(strcmp(Indicator,'NN') || ...
            strcmp(Indicator,'NN_Pwise') || ...
            strcmp(Indicator,'NN_modal_Pwise') || ...
            strcmp(Indicator,'NN_modal_patch_Pwise'))
        data_fname = sprintf('%s_IND_%s',data_fname,nn_model);
    else
        error('Indicator %s not available',Indicator);
    end
    if(~isempty(ind_var))
        data_fname = sprintf('%s_IVAR_%s',data_fname,ind_var);
    end
    data_fname = sprintf('%s_LIM_%s',data_fname,Limiter);
    if(~isempty(lim_var))
        data_fname = sprintf('%s_LVAR_%s',data_fname,lim_var);
    end
end