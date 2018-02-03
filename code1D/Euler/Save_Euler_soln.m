function Save_Euler_soln(q,gas_gamma, gas_const)

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
data_fname = sprintf('%s_LIM_%s.dat',data_fname,rec_limiter);


rho = q(:,:,1);
vel = q(:,:,2)./q(:,:,1);
pre = (gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
fid = fopen(data_fname,'w');
fprintf(fid, '%.16f %.16f %.16f  %.16f\n', [x(:) rho(:) vel(:) pre(:)]');
fclose(fid);


return