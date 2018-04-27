function fname = Euler_ind_fname1D(REL_PATH)

Globals1D_DG;
Globals1D_MLP;

fname = sprintf('%sUser_output/%s1D_%s_P%d_N%d',REL_PATH,model,test_name,...
                 N,Nelem);
if(strcmp(indicator_type,'minmod'))
    fname = sprintf('%s_IND_%s',fname,indicator_type);
elseif(strcmp(indicator_type,'TVB'))
    fname = sprintf('%s_IND_%s_%d',fname,indicator_type,TVB_M);
elseif(strcmp(indicator_type,'NN'))
    fname = sprintf('%s_IND_%s',fname,nn_model);
elseif(strcmp(indicator_type,'FuShu'))
    fname = sprintf('%s_IND_FS',fname);    
else
    error('Indicator %s not available',indicator_type);
end
fname = sprintf('%s_IVAR_%s_LIM_%s_LVAR_%s_tcells',...
                fname,ind_var,rec_limiter,lim_var);
if(mesh_pert ~= 0.0)
    fname = sprintf('%s_pert',fname);
end
fname = sprintf('%s.dat',fname); 


return
