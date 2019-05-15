function fname = Scalar_fnamebase1D(Problem,N,K,Limit,mesh_pert)

% Globals1D_DG;
% Globals1D_MLP;

mkdir('OUTPUT');

fname = sprintf('OUTPUT/%s1D_%s_P%d_N%d',Problem.model,Problem.test_name,N,K);
if(strcmp(Limit.Indicator,'NONE'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);
elseif(strcmp(Limit.Indicator,'ALL'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);    
elseif(strcmp(Limit.Indicator,'MINMOD'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);
elseif(strcmp(Limit.Indicator,'TVB'))
    fname = sprintf('%s_IND_%s_%d',fname,Limit.Indicator,Limit.TVBM);
elseif(strcmp(Limit.Indicator,'NN'))
    fname = sprintf('%s_IND_%s',fname,Limit.nn_model);
elseif(strcmp(Limit.Indicator,'FuShu'))
    fname = sprintf('%s_IND_FS',fname);    
else
    error('Indicator %s not available',Limit.Indicator);
end
fname = sprintf('%s_LIM_%s',fname,Limit.Limiter);
if(mesh_pert ~= 0.0)
    fname = sprintf('%s_pert',fname);
end


return
