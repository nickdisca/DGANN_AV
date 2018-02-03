function Save_SWE_soln(q)

Globals1D_DG;
Globals1D_MLP;

data_fname = sprintf('../User_output/%s1D_%s_P%d_N%d',model,test_name,N,...
                      Nelem);
if(strcmp(indicator_type,'minmod'))
    data_fname = sprintf('%s_IND_%s',data_fname,indicator_type);
elseif(strcmp(indicator_type,'TVB'))
    data_fname = sprintf('%s_IND_%s_%d',data_fname,indicator_type,TVB_M);
elseif(strcmp(indicator_type,'NN'))
    data_fname = sprintf('%s_IND_%s_%s_%s_%s_%s',data_fname,indicator_type,...
            nn_model,sub_model,data_set,data_subset);
else
    error('Indicator %s not available',indicator_type);
end
data_fname = sprintf('%s_IVAR_%s_LIM_%s_LVAR_%s.dat',...
                     data_fname,ind_var,rec_limiter,lim_var);


depth = q(:,:,1);
vel   = q(:,:,2)./q(:,:,1);
fid = fopen(data_fname,'w');
fprintf(fid, '%.16f %.16f %.16f \n', [x(:) depth(:) vel(:)]');
fclose(fid);


return