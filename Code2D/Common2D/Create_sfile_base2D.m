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

if(strcmp(Viscosity.model,'NONE'))
    data_fname = sprintf('%s_VISC_%s',data_fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'MDH'))
    %data_fname = sprintf('%s_VISC_%s_%d_%d_%d',fname,Viscosity.model,Viscosity.c_A,Viscosity.c_k,Viscosity.c_max);
    data_fname = sprintf('%s_VISC_%s',data_fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'MDA'))
    %data_fname = sprintf('%s_VISC_%s_%d',fname,Viscosity.model,Viscosity.c_max);
    data_fname = sprintf('%s_VISC_%s',data_fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'EV'))
    %data_fname = sprintf('%s_VISC_%s_%d_%d',fname,Viscosity.model,Viscosity.c_E,Viscosity.c_max);
    data_fname = sprintf('%s_VISC_%s',data_fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'NN'))
    %data_fname = sprintf('%s_VISC_%s',fname,Viscosity.nn_visc_model);
    data_fname = sprintf('%s_VISC_%s',data_fname,Viscosity.model);
else
    error('Viscosity model %s not available',Viscosity.model);
end