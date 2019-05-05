% Create save file base name (saved in th same folder)

data_fname = sprintf('%s2D_%s_P%d',Problem.model,Problem.test_name,Mesh.N);

if(strcmp(Limit.Limiter,'NONE'))
    data_fname = sprintf('%s_NONE',data_fname);
else
    if(strcmp(Limit.Indicator,'NONE') || strcmp(Limit.Indicator,'ALL'))
        data_fname = sprintf('%s_IND_%s',data_fname,Limit.Indicator);
    elseif(strcmp(Limit.Indicator,'TVB'))
        data_fname = sprintf('%s_IND_%s_%d',data_fname,Limit.Indicator,Limit.TVBM);   
    elseif(strcmp(Limit.Indicator,'NN'))
        data_fname = sprintf('%s_IND_%s',data_fname,Limit.nn_model);
    else
        error('Indicator %s not available',Limit.Indicator);
    end
    if(Limit.Filter_const)
        data_fname = sprintf('%s_ConstFilt',data_fname);
    end
%     if(Limit.Remove_iso)
%         data_fname = sprintf('%s_IsoFilt',data_fname);
%     end
    if(~isempty(Limit.ind_var))
        data_fname = sprintf('%s_IVAR_%s',data_fname,Limit.ind_var);
    end
    data_fname = sprintf('%s_LIM_%s',data_fname,Limit.Limiter);
    if(~isempty(Limit.lim_var))
        data_fname = sprintf('%s_LVAR_%s',data_fname,Limit.lim_var);
    end
end