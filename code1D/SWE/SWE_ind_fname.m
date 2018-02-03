function fname = SWE_ind_fname()

Globals1D_DG;
Globals1D_MLP;

fname = sprintf('../User_output/%s1D_%s_P%d_N%d',model,test_name,...
                 N,Nelem);
if(strcmp(indicator_type,'minmod'))
    fname = sprintf('%s_IND_%s',fname,indicator_type);
elseif(strcmp(indicator_type,'TVB'))
    fname = sprintf('%s_IND_%s_%d',fname,indicator_type,TVB_M);
elseif(strcmp(indicator_type,'NN'))
    fname = sprintf('%s_IND_%s_%s_%s_%s_%s',fname,indicator_type,...
            nn_model,sub_model,data_set,data_subset);
else
    error('Indicator %s not available',indicator_type);
end
fname = sprintf('%s_IVAR_%s_LIM_%s_LVAR_%s_tcells.dat',...
                fname,ind_var,rec_limiter,lim_var);


return