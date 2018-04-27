function Save_scalar_soln1D(u,REL_PATH)

Globals1D_DG;
Globals1D_MLP;

data_fname = sprintf('%sUser_output/%s1D_%s_P%d_N%d',REL_PATH,model,test_name,N,...
                      Nelem);
if(strcmp(indicator_type,'minmod'))
    data_fname = sprintf('%s_IND_%s',data_fname,indicator_type);
elseif(strcmp(indicator_type,'TVB'))
    data_fname = sprintf('%s_IND_%s_%d',data_fname,indicator_type,TVB_M);
elseif(strcmp(indicator_type,'NN'))
    data_fname = sprintf('%s_IND_%s',data_fname,nn_model);
elseif(strcmp(indicator_type,'FuShu'))
    data_fname = sprintf('%s_IND_FS',data_fname);
else    
    error('Indicator %s not available',indicator_type);
end
data_fname = sprintf('%s_LIM_%s',data_fname,rec_limiter);

if(mesh_pert ~= 0.0)
    data_fname = sprintf('%s_pert',data_fname);
end
data_fname = sprintf('%s.dat',data_fname);

fid = fopen(data_fname,'w');
fprintf(fid, '%.16f %.16f \n', [x(:) u(:)]');
fclose(fid);


return
